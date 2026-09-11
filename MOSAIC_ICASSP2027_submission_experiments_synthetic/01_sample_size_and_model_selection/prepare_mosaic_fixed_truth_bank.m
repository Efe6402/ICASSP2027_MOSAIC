function [truth,index]=prepare_mosaic_fixed_truth_bank(cfg)
%PREPARE_MOSAIC_FIXED_TRUTH_BANK Freeze the realized graph; vary signals only.
setup_three_method_paths();
if ~isfile(cfg.source_dataset), error('Source dataset not found:\n%s',cfg.source_dataset); end
if ~exist(cfg.data_dir,'dir'), mkdir(cfg.data_dir); end
truth_file=fullfile(cfg.data_dir,'MOSAIC_FIXED_TRUTH.mat');
if isfile(truth_file) && ~cfg.overwrite_data
    z=load(truth_file,'truth'); truth=z.truth;
else
    z=load(cfg.source_dataset);
    if isfield(z,'data'), d=z.data; elseif isfield(z,'truth'), d=z.truth;
    else, error('Source MAT must contain data or truth.'); end
    need={'n','K','r','A_state_true','B_true','Gamma_true','A_supra_true', ...
        'A_physical_true','A_copy_true','truth_masks', ...
        'truth_support_structural','Sigma','Sigma_clean','cfg'};
    for j=1:numel(need), assert(isfield(d,need{j}),'Missing source field %s.',need{j}); end
    truth=struct();
    for j=1:numel(need), truth.(need{j})=d.(need{j}); end
    optional={'Theta_true','edge_i','edge_j','W_contribution_realized', ...
        'copy_supports','community_nodes','topology_diagnostics','seed', ...
        'copy_scale_applied','realized_copy_inter_mass_fraction'};
    for j=1:numel(optional)
        if isfield(d,optional{j}), truth.(optional{j})=d.(optional{j}); end
    end
    truth.N=truth.n*truth.K;
    truth.A_within_true=truth.A_supra_true.*truth.truth_masks.within;
    truth.A_cross_true=truth.A_supra_true.*truth.truth_masks.crossnode_crossview;
    truth.A_copy_masked_true=truth.A_supra_true.*truth.truth_masks.copy;
    truth.source_dataset=cfg.source_dataset;
    save(truth_file,'truth','-v7.3');
    export_truth_artifacts(truth,cfg.data_dir);
end
verify_truth(truth);
index_file=fullfile(cfg.data_dir,'BANK_INDEX.csv');
if isfile(index_file)&&~cfg.overwrite_data
    index=readtable(index_file,'TextType','string');
    expected=2*numel(cfg.sample_counts)*cfg.validation_count;
    if cfg.validation_count~=cfg.test_count
        expected=numel(cfg.sample_counts)*(cfg.validation_count+cfg.test_count);
    end
    if height(index)==expected && all(ismember(cfg.sample_counts,unique(index.sample_count)))
        return
    end
end
index=generate_banks(truth,cfg);
writetable(index,index_file);
end

function index=generate_banks(truth,cfg)
split_col=strings(0,1); p_col=[]; trial_col=[]; seed_col=[]; file_col=strings(0,1);
splits=["validation","test"]; counts=[cfg.validation_count,cfg.test_count];
[V,D]=eig(0.5*(truth.Sigma+truth.Sigma'));
ev=max(real(diag(D)),0); F=real(V*diag(sqrt(ev)));
maxp=max(cfg.sample_counts);
for si=1:2
    for trial=1:counts(si)
        seed=cfg.seed_base+si*100000+trial;
        old=rng; clean=onCleanup(@()rng(old)); rng(seed,'twister');
        Ymax=F*randn(truth.N,maxp);
        for pi=1:numel(cfg.sample_counts)
            p=cfg.sample_counts(pi); Y=Ymax(:,1:p);
            record=make_record(Y,truth,p,char(splits(si)),trial,seed);
            out=fullfile(cfg.data_dir,sprintf('%s_p%05d_trial%02d.mat', ...
                splits(si),p,trial));
            save(out,'record','-v7.3');
            split_col(end+1,1)=splits(si); p_col(end+1,1)=p; %#ok<AGROW>
            trial_col(end+1,1)=trial; seed_col(end+1,1)=seed; %#ok<AGROW>
            file_col(end+1,1)=string(out); %#ok<AGROW>
        end
        clear clean Ymax
    end
end
index=table(split_col,p_col,trial_col,seed_col,file_col, ...
    'VariableNames',{'split','sample_count','trial','seed','file'});
end

function record=make_record(Y,truth,p,split,trial,seed)
Cy=(Y*Y')/p; d=diag(Cy); S=max(d+d'-2*Cy,0);
S=0.5*(S+S'); S(1:size(S,1)+1:end)=0;
Sp=zeros(truth.n); Sq=zeros(truth.K);
for s=1:p
    Xi=reshape(Y(:,s),truth.n,truth.K);
    Sp=Sp+cov(Xi'); Sq=Sq+cov(Xi);
end
record=struct('split',split,'sample_count',p,'trial',trial,'seed',seed, ...
    'n',truth.n,'K',truth.K,'S',S,'Cy_uncentered',Cy, ...
    'Sp',Sp/p,'Sq',Sq/p);
end

function export_truth_artifacts(t,out)
writematrix(t.A_state_true(:,:,1),fullfile(out,'A_MODE_1_BINARY.csv'));
writematrix(t.A_state_true(:,:,2),fullfile(out,'A_MODE_2_BINARY.csv'));
writematrix(t.B_true,fullfile(out,'B_TRUE_COLUMN_SIMPLEX.csv'));
writematrix(t.Gamma_true,fullfile(out,'GAMMA_TRUE.csv'));
writematrix(t.A_supra_true,fullfile(out,'A_SUPRA_TRUE_CONTINUOUS.csv'));
writematrix(t.A_physical_true,fullfile(out,'A_PHYSICAL_TRUE_CONTINUOUS.csv'));
writematrix(t.A_copy_true,fullfile(out,'A_COPY_TRUE_CONTINUOUS.csv'));
writematrix(t.truth_support_structural.all,fullfile(out,'TRUTH_SUPPORT_FULL_BINARY.csv'));
writematrix(t.truth_support_structural.within,fullfile(out,'TRUTH_SUPPORT_WITHIN_BINARY.csv'));
writematrix(t.truth_support_structural.crossnode_crossview,fullfile(out,'TRUTH_SUPPORT_CROSS_BINARY.csv'));
writematrix(t.truth_support_structural.copy,fullfile(out,'TRUTH_SUPPORT_COPY_BINARY.csv'));
writematrix(t.Sigma,fullfile(out,'SIGNAL_COVARIANCE_SIGMA.csv'));
if isfield(t,'W_contribution_realized')
    writematrix(t.W_contribution_realized,fullfile(out,'VIEW_MODAL_CONTRIBUTIONS.csv'));
end
end

function verify_truth(t)
assert(t.N==t.n*t.K&&isequal(size(t.A_supra_true),[t.N,t.N]));
assert(norm(t.A_supra_true-t.A_supra_true','fro')<1e-10);
assert(max(abs(diag(t.A_supra_true)))<1e-10);
assert(all(t.B_true(:)>=0)&&max(abs(sum(t.B_true,1)-1))<1e-10);
assert(isequal(t.A_supra_true>0,t.truth_support_structural.all));
end
