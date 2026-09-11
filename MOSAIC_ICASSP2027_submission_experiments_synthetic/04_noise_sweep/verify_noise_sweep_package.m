function report=verify_noise_sweep_package(cfg)
%VERIFY_NOISE_SWEEP_PACKAGE Check self-contained files and frozen selections.
if nargin<1||isempty(cfg), cfg=noise_sweep_config('quick'); end
required={'RUN_FULL_NOISE_SWEEP.m','RUN_QUICK_TEST.m','noise_sweep_config.m', ...
    'view_variance_profile_from_kappa.m','prepare_noise_signal_bank.m', ...
    'run_noise_sweep.m','fit_fixed_candidate.m','plot_noise_results.m', ...
    'graph_auc_metrics.m','apply_graph_f1_thresholds.m','make_signal_record.m', ...
    'load_default_truth.m','setup_noise_sweep_paths.m','validate_noise_bank.m', ...
    'signal_covariance_from_supra.m'};
for j=1:numel(required)
    assert(isfile(fullfile(cfg.root,required{j})),'Missing file: %s',required{j});
end
assert(isfile(cfg.default_truth_file),'Missing source truth MAT file.');
expected=struct('MOSAIC',[407 420 421], ...
    'PGL2021',[269 294 319],'ZW2024',[593 482 592]);
methods=string(fieldnames(expected));
for method=methods'
    file=fullfile(cfg.selected_config_dir,method+"_F1_COMPONENT_TOP3.csv");
    assert(isfile(file),'Missing selected configuration file: %s',file);
    T=readtable(file,'TextType','string');
    T=T(T.selection_criterion=="f1_component",:); T=sortrows(T,'selection_rank');
    assert(isequal(double(T.candidate_id(:))',expected.(char(method))));
    assert(all(isfinite(double(T{:,{'tau_supra','tau_within','tau_cross','tau_copy'}})),'all'));
end
solverfiles={fullfile(cfg.root,'solver','mosaic','solve_mosaic_crossview.m'), ...
    fullfile(cfg.root,'solver','pgl2021','Learn_PGL.m'), ...
    fullfile(cfg.root,'solver','zw2024','solve_zhang_wai_2024_paper.m')};
assert(all(cellfun(@isfile,solverfiles)),'At least one package-local solver is missing.');
truth=load_default_truth(cfg);
assert(truth.n==16&&truth.K==4&&truth.r==2&&truth.N==64);
assert(max(abs(sum(truth.B_true,1)-1))<1e-12);
if strcmpi(cfg.mode,'full')
    assert(isequal(cfg.snr_db,[Inf 20 15 12 10 8 6 4 2 0]));
    assert(isequal(cfg.kappa_grid,[0 0.2 0.4 0.6 0.8 1.0]));
end
assert(isequal(cfg.view_variance_endpoint,[0.4 0.8 1.2 1.6]));
assert(abs(mean(cfg.view_variance_endpoint)-1)<1e-12);
assert(isequal(view_variance_profile_from_kappa(0,cfg.view_variance_endpoint),ones(1,4)));
assert(max(abs(view_variance_profile_from_kappa(1,cfg.view_variance_endpoint)- ...
    cfg.view_variance_endpoint))<1e-12);
report=struct('valid',true,'method_count',3,'top_k_available',3, ...
    'source_p',5000,'n',truth.n,'K',truth.K,'r',truth.r, ...
    'noise_model','viewwise_kappa_heteroscedastic_gaussian', ...
    'view_dependent_only',true,'node_dependent_noise',false, ...
    'kappa_zero_is_equal_variance_iid',true, ...
    'kappa_one_is_original_view_variance_endpoint',true);
end
