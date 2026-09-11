function summary = run_coil20_experiment_core(package_root)
%RUN_COIL20_EXPERIMENT_CORE Six-view COIL-20 selection and final fit.

if nargin < 1 || isempty(package_root)
    package_root = fileparts(mfilename('fullpath'));
end

cfg = experiment_configuration(package_root);
addpath(fullfile(package_root,'solver'));

selection_dir = fullfile(cfg.output_dir,'selection');
final_dir = fullfile(cfg.output_dir,'final');
paper_dir = fullfile(package_root,'paper_figures');
ensure_dir(selection_dir); ensure_dir(final_dir); ensure_dir(paper_dir);

[target_data,train_data,validation_data] = load_experiment_data(cfg);
profiles = focused_profile_bank();
writetable(struct2table(profiles,'AsArray',true), ...
    fullfile(selection_dir,'hyperparameter_bank.csv'));

fprintf('\nCOIL-20 six-view auxiliary cross-configuration selection\n');
fprintf('Ranks: %s | profiles: %d | source/validation pairs: %d\n', ...
    mat2str(cfg.rank_candidates),numel(profiles),numel(train_data));

[aggregate_table,rank_table] = run_selection( ...
    train_data,validation_data,profiles,cfg);
writetable(aggregate_table,fullfile(selection_dir, ...
    'profile_rank_validation.csv'));

[selected_rank,knee_status,knee_score] = select_rank(rank_table);
rank_table.knee_score = knee_score;
rank_table.selected = rank_table.r == selected_rank;
writetable(rank_table,fullfile(selection_dir,'rank_selection.csv'));
save_rank_curve(rank_table,fullfile(selection_dir,'rank_validation_curve.png'));

at_rank = aggregate_table(aggregate_table.r == selected_rank & ...
    aggregate_table.valid,:);
at_rank = sortrows(at_rank,{'mean_validation_score','profile_id'}, ...
    {'ascend','ascend'});
if isempty(at_rank)
    error('No valid profile is available at selected rank %d.',selected_rank);
end
selected_row = at_rank(1,:);
selected_profile = profiles([profiles.profile_id] == selected_row.profile_id);
selected_table = selected_row;
selected_table.selection_rank = 1;
selected_table.knee_status = string(knee_status);
selected_table.final_seed = cfg.seed + cfg.final_seed_offset + ...
    cfg.case_seed + 10000*selected_rank;
writetable(selected_table,fullfile(selection_dir,'SELECTED_TOP1.csv'));

fprintf('\nSelected r=%d, profile %d (%s), auxiliary score %.15g.\n', ...
    selected_rank,selected_profile.profile_id,selected_profile.name, ...
    selected_row.mean_validation_score);

[compactFit,analysis,final_record] = fit_final_target( ...
    target_data,selected_profile,selected_rank,selected_table.final_seed,cfg);
candidate_id = sprintf('top01_r%d_profile%d_%s',selected_rank, ...
    selected_profile.profile_id,selected_profile.name);
fit_file = fullfile(final_dir,[candidate_id '.mat']);
profile = selected_profile; %#ok<NASGU>
record = final_record; %#ok<NASGU>
save(fit_file,'compactFit','analysis','profile','record','-v7.3');

final_table = struct2table(final_record,'AsArray',true);
final_table.auxiliary_validation_score = ...
    repmat(selected_row.mean_validation_score,height(final_table),1);
final_table.modal_rank = repmat(selected_rank,height(final_table),1);
final_table.profile_id = repmat(selected_profile.profile_id,height(final_table),1);
final_table.profile_name = repmat(string(selected_profile.name),height(final_table),1);
writetable(final_table,fullfile(final_dir,'FINAL_TOP1_FIT.csv'));
writematrix(compactFit.B,fullfile(final_dir,'B_modal_view_profiles.csv'));
writematrix(compactFit.view_graph,fullfile(final_dir,'induced_view_graph.csv'));

generate_coil20_paper_figures(fit_file,cfg.target_file,paper_dir);

expected_rank = 3; expected_profile = 2508;
matches_packaged = selected_rank == expected_rank && ...
    selected_profile.profile_id == expected_profile;
if ~matches_packaged
    warning(['Selected configuration differs from the packaged paper result. ' ...
        'Inspect results/selection before using the final figures.']);
end

summary = struct();
summary.output_dir = cfg.output_dir;
summary.selected_rank = selected_rank;
summary.selected_profile_id = selected_profile.profile_id;
summary.selected_profile_name = selected_profile.name;
summary.mean_auxiliary_validation_score = selected_row.mean_validation_score;
summary.matches_packaged_top1 = matches_packaged;
summary.final_fit_file = fit_file;
summary.paper_figure_dir = paper_dir;
save(fullfile(cfg.output_dir,'RUN_SUMMARY.mat'),'summary','cfg','-v7.3');
end

function cfg = experiment_configuration(package_root)
cfg = struct();
cfg.seed = 20260811;
cfg.case_seed = 110000;
cfg.final_seed_offset = 900000;
cfg.rank_candidates = 2:6;
cfg.n_restarts_selection = 2;
cfg.n_restarts_final = 2;
cfg.max_iter_selection = 600;
cfg.max_iter_final = 1200;
cfg.stationarity_tol = 1e-6;
cfg.objective_tol = 1e-9;
cfg.output_dir = fullfile(package_root,'results');

data_root = fullfile(package_root,'data','global_6view_phases');
cfg.target_file = fullfile(data_root,'coil20_global6_phase_00.mat');
cfg.train_files = { ...
    fullfile(data_root,'coil20_global6_phase_05.mat'), ...
    fullfile(data_root,'coil20_global6_phase_25.mat'), ...
    fullfile(data_root,'coil20_global6_phase_45.mat')};
cfg.validation_files = { ...
    fullfile(data_root,'coil20_global6_phase_10.mat'), ...
    fullfile(data_root,'coil20_global6_phase_30.mat'), ...
    fullfile(data_root,'coil20_global6_phase_50.mat')};
end

function [target_data,train_data,validation_data] = load_experiment_data(cfg)
all_files = [{cfg.target_file},cfg.train_files,cfg.validation_files];
for q = 1:numel(all_files)
    if ~isfile(all_files{q}), error('Missing data file: %s',all_files{q}); end
end

target_data = load_one(cfg.target_file);
train_data = cellfun(@load_one,cfg.train_files,'UniformOutput',false);
validation_data = cellfun(@load_one,cfg.validation_files,'UniformOutput',false);
for q = 1:numel(train_data)
    assert_aligned(target_data,train_data{q},cfg.train_files{q});
    assert_aligned(target_data,validation_data{q},cfg.validation_files{q});
end
end

function data = load_one(file_path)
S = load(file_path,'data');
if ~isfield(S,'data'), error('Variable data is absent from %s.',file_path); end
data = S.data;
required = {'X','view_angles_deg','object_ids'};
for q = 1:numel(required)
    if ~isfield(data,required{q}), error('data.%s is required.',required{q}); end
end
data.num_views = numel(data.X);
data.num_objects = size(data.X{1},1);
data.signal_dimension = size(data.X{1},2);
if data.num_views ~= 6 || data.num_objects ~= 20 || data.signal_dimension ~= 1024
    error('Expected a 20-object, six-view, 1024-dimensional COIL dataset.');
end
for k = 1:data.num_views
    if ~isequal(size(data.X{k}),[data.num_objects,data.signal_dimension]) || ...
            any(~isfinite(data.X{k}(:)))
        error('Invalid signal matrix in view %d.',k);
    end
end
data.object_ids = double(data.object_ids(:));
data.view_angles_deg = double(data.view_angles_deg(:)');
end

function assert_aligned(reference,data,label)
if data.num_views ~= reference.num_views || ...
        data.num_objects ~= reference.num_objects || ...
        data.signal_dimension ~= reference.signal_dimension || ...
        ~isequal(data.object_ids,reference.object_ids)
    error('Dataset %s is not aligned with the phase-0 target.',label);
end
end

function [aggregate_table,rank_table] = run_selection(train_data,val_data,profiles,cfg)
n_sources = numel(train_data);
n_profiles = numel(profiles);
n_ranks = numel(cfg.rank_candidates);

Ytrain = cell(n_sources,1); Strain = cell(n_sources,1); Sval = cell(n_sources,1);
for s = 1:n_sources
    Ytrain{s} = vertcat(train_data{s}.X{:});
    Strain{s} = distance_matrix(Ytrain{s});
    Yval = vertcat(val_data{s}.X{:});
    Sval{s} = distance_matrix(Yval);
end

template = struct('r',NaN,'profile_id',NaN,'profile_name','', ...
    'valid',false,'n_validation_scores',0,'mean_validation_score',NaN, ...
    'std_validation_score',NaN,'se_validation_score',NaN, ...
    'min_validation_score',NaN,'max_validation_score',NaN, ...
    'mean_within_score',NaN,'mean_cross_score',NaN,'mean_copy_score',NaN);
rows = repmat(template,n_ranks*n_profiles,1);
row_index = 0;

for ir = 1:n_ranks
    r = cfg.rank_candidates(ir);
    for ip = 1:n_profiles
        p = profiles(ip);
        scores = nan(n_sources,1); within = scores; cross = scores; copy = scores;
        fprintf('[selection] r=%d | profile %d/%d: %s\n', ...
            r,ip,n_profiles,p.name);

        for s = 1:n_sources
            paired_seed = cfg.seed + cfg.case_seed + 1000000*s + 10000*r;
            observed = struct('n',train_data{s}.num_objects, ...
                'K',train_data{s}.num_views,'r',r,'Y',Ytrain{s}, ...
                'S',Strain{s},'cost_scale',1);
            opts = profile_options(p,paired_seed,cfg,'selection');
            try
                fit = solve_mosaic_crossview(observed,opts);
                vm = validation_metrics(fit.A_supra,Sval{s}, ...
                    train_data{s}.num_objects,train_data{s}.num_views);
                scores(s) = vm.total_score;
                within(s) = vm.within_score;
                cross(s) = vm.cross_score;
                copy(s) = vm.copy_score;
            catch ME
                warning('Selection fit failed for source %d, r=%d, profile %d: %s', ...
                    s,r,p.profile_id,ME.message);
            end
        end

        row_index = row_index + 1;
        row = template;
        row.r = r; row.profile_id = p.profile_id; row.profile_name = p.name;
        good = isfinite(scores);
        row.n_validation_scores = nnz(good);
        row.valid = all(good);
        if row.valid
            row.mean_validation_score = mean(scores);
            row.std_validation_score = std(scores,0);
            row.se_validation_score = row.std_validation_score/sqrt(n_sources);
            row.min_validation_score = min(scores);
            row.max_validation_score = max(scores);
            row.mean_within_score = mean(within);
            row.mean_cross_score = mean(cross);
            row.mean_copy_score = mean(copy);
        end
        rows(row_index) = row;
    end
end

aggregate_table = struct2table(rows,'AsArray',true);
aggregate_table = sortrows(aggregate_table,{'r','mean_validation_score'}, ...
    {'ascend','ascend'});

rank_template = struct('r',NaN,'best_profile_id',NaN, ...
    'best_profile_name','','validation_score',NaN,'validation_std',NaN, ...
    'validation_se',NaN,'incremental_gain',NaN,'relative_gain',NaN, ...
    'knee_score',NaN);
rank_rows = repmat(rank_template,n_ranks,1);
for ir = 1:n_ranks
    r = cfg.rank_candidates(ir);
    sub = aggregate_table(aggregate_table.r == r & aggregate_table.valid,:);
    sub = sortrows(sub,{'mean_validation_score','profile_id'},{'ascend','ascend'});
    if isempty(sub), error('No complete valid candidate exists for r=%d.',r); end
    rank_rows(ir).r = r;
    rank_rows(ir).best_profile_id = sub.profile_id(1);
    rank_rows(ir).best_profile_name = sub.profile_name{1};
    rank_rows(ir).validation_score = sub.mean_validation_score(1);
    rank_rows(ir).validation_std = sub.std_validation_score(1);
    rank_rows(ir).validation_se = sub.se_validation_score(1);
    if ir > 1
        rank_rows(ir).incremental_gain = rank_rows(ir-1).validation_score - ...
            rank_rows(ir).validation_score;
        rank_rows(ir).relative_gain = rank_rows(ir).incremental_gain / ...
            max(abs(rank_rows(ir-1).validation_score),eps);
    end
end
rank_table = struct2table(rank_rows,'AsArray',true);
end

function [selected_rank,status,knee] = select_rank(rank_table)
r = double(rank_table.r(:)); V = double(rank_table.validation_score(:));
x = (r-r(1))/max(r(end)-r(1),eps);
total_drop = V(1)-V(end);
if ~isfinite(total_drop) || total_drop <= 1e-10*max(1,abs(V(1)))
    knee = zeros(size(r)); [~,index] = min(V);
    status = 'fallback_minimum_no_overall_decrease';
else
    y = (V(1)-V)/total_drop;
    knee = y-x; knee(1)=0; knee(end)=0;
    interior = 2:numel(r)-1;
    [maximum,j] = max(knee(interior)); index = interior(j);
    if ~isfinite(maximum) || maximum <= 0
        [~,index] = min(V); status = 'fallback_minimum_no_positive_knee';
    else
        status = 'geometric_endpoint_chord_knee';
    end
end
selected_rank = r(index);
end

function [compactFit,analysis,record] = fit_final_target(data,p,r,seed,cfg)
Y = vertcat(data.X{:});
observed = struct('n',data.num_objects,'K',data.num_views,'r',r, ...
    'Y',Y,'S',distance_matrix(Y),'cost_scale',1);
opts = profile_options(p,seed,cfg,'final');
timer = tic;
fit = solve_mosaic_crossview(observed,opts);
runtime_seconds = toc(timer);
compactFit = compact_fit(fit);
analysis = analyze_fit(compactFit,data);
record = struct('valid',true,'converged',strcmpi(fit.status,'converged'), ...
    'status',fit.status,'seed',seed,'runtime_seconds',runtime_seconds, ...
    'objective',fit.objective,'best_restart',fit.best_restart, ...
    'stationarity',fit.stationarity_combined,'iterations',fit.iterations);
end

function opts = profile_options(p,seed,cfg,stage)
opts = struct('seed',seed,'stationarity_tol',cfg.stationarity_tol, ...
    'objective_tol',cfg.objective_tol,'verbose',false, ...
    'alpha_a',p.alpha_a,'alpha_c',p.alpha_c, ...
    'beta_a',p.beta_a,'beta_c',p.beta_c, ...
    'tau_intra',p.tau_intra,'tau_inter',p.tau_inter, ...
    'eta_c',p.eta_c,'lambda_B',p.lambda_B, ...
    'epsilon_intra',p.epsilon_intra,'epsilon_inter',p.epsilon_inter);
if strcmp(stage,'selection')
    opts.n_restarts = cfg.n_restarts_selection;
    opts.max_iter = cfg.max_iter_selection;
else
    opts.n_restarts = cfg.n_restarts_final;
    opts.max_iter = cfg.max_iter_final;
end
end

function S = distance_matrix(Y)
p = size(Y,2);
d = sum(Y.^2,2)/p;
S = max(d+d'-2*(Y*Y'/p),0);
S = 0.5*(S+S'); S(1:size(S,1)+1:end)=0;
values = S(triu(true(size(S)),1));
values = values(isfinite(values) & values>0);
if isempty(values), error('All positive pairwise distances vanished.'); end
S = S/median(values);
end

function vm = validation_metrics(A,S,n,K)
N = n*K; A = max(0,0.5*(A+A')); A(1:N+1:end)=0;
upper = triu(true(N),1); mass = sum(A(upper));
if ~isfinite(mass) || mass<=0
    vm = struct('total_score',Inf,'within_score',Inf, ...
        'cross_score',Inf,'copy_score',Inf); return;
end
A = A/mass;
within=false(N); cross=false(N); copy=false(N);
for k=1:K
    Ik=(k-1)*n+(1:n); within(Ik,Ik)=~eye(n);
    for ell=k+1:K
        Il=(ell-1)*n+(1:n);
        cross_block=true(n); cross_block(1:n+1:end)=false;
        cross(Ik,Il)=cross_block;
        copy_block=false(n); copy_block(1:n+1:end)=true;
        copy(Ik,Il)=copy_block;
    end
end
within=within&upper; cross=cross&upper; copy=copy&upper;
vm = struct('total_score',sum(A(upper).*S(upper)), ...
    'within_score',sum(A(within).*S(within)), ...
    'cross_score',sum(A(cross).*S(cross)), ...
    'copy_score',sum(A(copy).*S(copy)));
end

function compact = compact_fit(fit)
compact = fit;
if isfield(compact,'observed'), compact=rmfield(compact,'observed'); end
if isfield(compact,'all_restarts')
    rr=compact.all_restarts;
    summary=repmat(struct('restart',NaN,'status','','objective',NaN, ...
        'stationarity_combined',NaN,'iterations',NaN),numel(rr),1);
    for q=1:numel(rr)
        summary(q).restart=q; summary(q).status=rr{q}.status;
        summary(q).objective=rr{q}.objective;
        summary(q).stationarity_combined=rr{q}.stationarity_combined;
        summary(q).iterations=rr{q}.iterations;
    end
    compact.restart_summary=summary;
    compact=rmfield(compact,'all_restarts');
end
end

function analysis = analyze_fit(fit,data)
n=data.num_objects; K=data.num_views; r=size(fit.B,2);
within=zeros(n,n,K);
for k=1:K
    for m=1:r
        within(:,:,k)=within(:,:,k)+fit.B(k,m)^2*fit.A_state(:,:,m);
    end
end
analysis = struct('within_view_adjacencies',within);
end

function save_rank_curve(T,file_path)
f=figure('Visible','off','Color','w','Position',[100 100 1050 650]);
plot(T.r,T.validation_score,'-o','LineWidth',2,'MarkerSize',8); hold on;
grid on; box on; xlabel('modal rank r');
ylabel('minimum mean held-out smoothness');
title('COIL-20 six-view auxiliary rank selection');
exportgraphics(f,file_path,'Resolution',240); close(f);
end

function profiles = focused_profile_bank()
rows = { ...
2501,'G01_anchor_profile25_exact',.003,.0003,.08,.020,.08,.30,.50,.0050,1e-5; ...
2502,'G02_eta_c_030',.003,.0003,.08,.020,.08,.30,.30,.0050,1e-5; ...
2503,'G03_eta_c_040',.003,.0003,.08,.020,.08,.30,.40,.0050,1e-5; ...
2504,'G04_eta_c_060',.003,.0003,.08,.020,.08,.30,.60,.0050,1e-5; ...
2505,'G05_eta_c_075',.003,.0003,.08,.020,.08,.30,.75,.0050,1e-5; ...
2506,'G06_eta_c_090',.003,.0003,.08,.020,.08,.30,.90,.0050,1e-5; ...
2507,'G07_eta_c_120',.003,.0003,.08,.020,.08,.30,1.20,.0050,1e-5; ...
2508,'G08_beta_c_010',.003,.0003,.08,.010,.08,.30,.50,.0050,1e-5; ...
2509,'G09_beta_c_015',.003,.0003,.08,.015,.08,.30,.50,.0050,1e-5; ...
2510,'G10_beta_c_025',.003,.0003,.08,.025,.08,.30,.50,.0050,1e-5; ...
2511,'G11_beta_c_030',.003,.0003,.08,.030,.08,.30,.50,.0050,1e-5; ...
2512,'G12_beta_c_040',.003,.0003,.08,.040,.08,.30,.50,.0050,1e-5; ...
2513,'G13_tau_inter_020',.003,.0003,.08,.020,.08,.20,.50,.0050,1e-5; ...
2514,'G14_tau_inter_025',.003,.0003,.08,.020,.08,.25,.50,.0050,1e-5; ...
2515,'G15_tau_inter_035',.003,.0003,.08,.020,.08,.35,.50,.0050,1e-5; ...
2516,'G16_tau_inter_040',.003,.0003,.08,.020,.08,.40,.50,.0050,1e-5; ...
2517,'G17_alpha_c_0001',.003,.0001,.08,.020,.08,.30,.50,.0050,1e-5; ...
2518,'G18_alpha_c_0002',.003,.0002,.08,.020,.08,.30,.50,.0050,1e-5; ...
2519,'G19_alpha_c_0005',.003,.0005,.08,.020,.08,.30,.50,.0050,1e-5; ...
2520,'G20_alpha_c_0008',.003,.0008,.08,.020,.08,.30,.50,.0050,1e-5; ...
2521,'G21_lambda_B_0000',.003,.0003,.08,.020,.08,.30,.50,.0000,1e-5; ...
2522,'G22_lambda_B_0025',.003,.0003,.08,.020,.08,.30,.50,.0025,1e-5; ...
2523,'G23_lambda_B_0100',.003,.0003,.08,.020,.08,.30,.50,.0100,1e-5; ...
2524,'G24_combo_eta075_beta015',.003,.0003,.08,.015,.08,.30,.75,.0050,1e-5; ...
2525,'G25_combo_eta060_beta015_tau035',.003,.0003,.08,.015,.08,.35,.60,.0050,1e-5};
template=struct('profile_id',NaN,'name','','alpha_a',NaN,'alpha_c',NaN, ...
    'beta_a',NaN,'beta_c',NaN,'tau_intra',NaN,'tau_inter',NaN, ...
    'eta_c',NaN,'lambda_B',NaN,'epsilon_intra',NaN,'epsilon_inter',NaN);
profiles=repmat(template,size(rows,1),1);
for q=1:size(rows,1)
    profiles(q).profile_id=rows{q,1}; profiles(q).name=rows{q,2};
    profiles(q).alpha_a=rows{q,3}; profiles(q).alpha_c=rows{q,4};
    profiles(q).beta_a=rows{q,5}; profiles(q).beta_c=rows{q,6};
    profiles(q).tau_intra=rows{q,7}; profiles(q).tau_inter=rows{q,8};
    profiles(q).eta_c=rows{q,9}; profiles(q).lambda_B=rows{q,10};
    profiles(q).epsilon_intra=rows{q,11}; profiles(q).epsilon_inter=rows{q,11};
end
end

function ensure_dir(folder)
if ~isfolder(folder), mkdir(folder); end
end
