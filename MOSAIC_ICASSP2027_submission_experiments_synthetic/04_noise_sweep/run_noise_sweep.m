function summary=run_noise_sweep(cfg)
%RUN_NOISE_SWEEP Evaluate frozen models over view-wise kappa x SNR noise.
setup_noise_sweep_paths(cfg.root);
verify_noise_sweep_package(cfg);
if isfolder(cfg.output_dir), error('Output already exists: %s',cfg.output_dir); end
mkdir(cfg.output_dir);
[truth,index,levels]=prepare_noise_signal_bank(cfg);
bank_validation=validate_noise_bank(cfg,truth,index,levels);
writetable(levels,fullfile(cfg.output_dir,'NOISE_CONDITIONS.csv'));
writetable(index,fullfile(cfg.output_dir,'BANK_INDEX.csv'));
writetable(struct2table(bank_validation),fullfile(cfg.output_dir,'BANK_VALIDATION.csv'));
save(fullfile(cfg.output_dir,'RUN_CONFIG.mat'),'cfg','-v7.3');
diary(fullfile(cfg.output_dir,'FULL_CONSOLE_LOG.txt'));
cleaner=onCleanup(@()diary('off'));

methods=["MOSAIC","PGL2021","ZW2024"];
allsummary=table(); allreal=table(); tauaudit=table(); fitstatus=table();
for method=methods
    selected=load_selected(cfg,method);
    writetable(selected,fullfile(cfg.output_dir,method+"_FROZEN_CONFIGURATIONS.csv"));
    for q=1:height(selected)
        tau=tau_from_row(selected(q,:));
        tauaudit=append_rows(tauaudit,tau_row(method,selected(q,:),tau));
        for lev=1:height(levels)
            Q=index(index.level==lev,:);
            X=nan(height(Q),10); runtime=nan(height(Q),1); status=strings(height(Q),1);
            for j=1:height(Q)
                [est,fitfile]=load_or_fit(cfg,method,selected(q,:),Q(j,:),truth,levels(lev,:));
                a=graph_auc_metrics(est.A_primary,truth);
                f=apply_graph_f1_thresholds(est.A_primary,truth,tau);
                X(j,:)=[a.supra,a.within,a.cross,a.copy,a.component, ...
                    f.supra,f.within,f.cross,f.copy,f.component];
                runtime(j)=est.runtime; status(j)=string(est.status);
                row=base_row(cfg,method,selected(q,:),levels(lev,:),Q.trial(j),tau);
                row=add_realization_metrics(row,X(j,:));
                row.realized_snr_db=Q.realized_snr_db(j);
                row.realized_noise_to_signal_power_ratio=Q.realized_noise_to_signal_power_ratio(j);
                row.realized_noise_energy_fraction=Q.realized_noise_energy_fraction(j);
                row.variance_profile_rotation=Q.variance_profile_rotation(j);
                row.applied_view_variance_1=Q.applied_view_variance_1(j);
                row.applied_view_variance_2=Q.applied_view_variance_2(j);
                row.applied_view_variance_3=Q.applied_view_variance_3(j);
                row.applied_view_variance_4=Q.applied_view_variance_4(j);
                row.realized_view_noise_power_1=Q.realized_view_noise_power_1(j);
                row.realized_view_noise_power_2=Q.realized_view_noise_power_2(j);
                row.realized_view_noise_power_3=Q.realized_view_noise_power_3(j);
                row.realized_view_noise_power_4=Q.realized_view_noise_power_4(j);
                row.observed_power=Q.observed_power(j);
                row.runtime_seconds=runtime(j); row.status=status(j);
                allreal=append_rows(allreal,row);
                fitstatus=append_rows(fitstatus,table(method,double(selected.candidate_id(q)), ...
                    double(levels.kappa(lev)),double(levels.target_snr_db(lev)), ...
                    lev,Q.trial(j),status(j),runtime(j),string(fitfile), ...
                    'VariableNames',{'method','candidate_id','kappa','target_snr_db', ...
                    'level','trial','status','runtime_seconds','fit_file'}));
                writetable(fitstatus,fullfile(cfg.output_dir,'FIT_STATUS_CHECKPOINT.csv'));
                if cfg.verbose
                    fprintf('%s candidate %d | kappa %.2f | SNR %s dB | trial %d | %s\n', ...
                        method,selected.candidate_id(q),levels.kappa(lev), ...
                        snr_label(levels.target_snr_db(lev)),Q.trial(j),status(j));
                end
            end
            row=base_row(cfg,method,selected(q,:),levels(lev,:),NaN,tau);
            row=add_summary_metrics(row,X);
            row.mean_runtime_seconds=mean(runtime,'omitnan');
            row.converged_fraction=mean(status=="converged"|status=="completed");
            allsummary=append_rows(allsummary,row);
            writetable(allsummary,fullfile(cfg.output_dir,'TEST_SUMMARY_CHECKPOINT.csv'));
            writetable(allreal,fullfile(cfg.output_dir,'TEST_REALIZATIONS_CHECKPOINT.csv'));
        end
    end
end
writetable(fitstatus,fullfile(cfg.output_dir,'FIT_STATUS.csv'));
writetable(allsummary,fullfile(cfg.output_dir,'TEST_SUMMARY.csv'));
writetable(allreal,fullfile(cfg.output_dir,'TEST_REALIZATIONS.csv'));
writetable(tauaudit,fullfile(cfg.output_dir,'FROZEN_THRESHOLD_AUDIT.csv'));
if cfg.make_figures, plot_noise_results(allsummary,cfg,cfg.output_dir); end
save(fullfile(cfg.output_dir,'COMPLETE_RESULTS.mat'),'cfg','levels','bank_validation', ...
    'allsummary','allreal','tauaudit','fitstatus','-v7.3');
write_protocol(cfg,levels);
summary=struct('output_dir',cfg.output_dir,'noise_conditions',levels, ...
    'test_summary',allsummary,'test_realizations',allreal, ...
    'frozen_thresholds',tauaudit,'fit_status',fitstatus);
clear cleaner
end

function C=load_selected(cfg,method)
file=fullfile(cfg.selected_config_dir,method+"_F1_COMPONENT_TOP3.csv");
assert(isfile(file),'Missing frozen configuration table: %s',file);
C=readtable(file,'TextType','string');
C=C(C.selection_criterion=="f1_component"&C.selection_rank<=cfg.top_k,:);
C=sortrows(C,'selection_rank'); assert(height(C)==cfg.top_k);
end

function [est,file]=load_or_fit(cfg,method,row,bankrow,truth,L)
% Clean data are identical for all kappa values within a trial, so cache the
% clean fit once.  Finite-SNR conditions are cached separately by kappa/SNR.
if isinf(double(L.target_snr_db))
    condition_folder='clean_shared_across_kappa';
else
    condition_folder=sprintf('kappa_%03d_snr_%02d', ...
        round(100*double(L.kappa)),round(double(L.target_snr_db)));
end
folder=fullfile(cfg.output_dir,'FIT_CACHE',char(method),condition_folder);
if ~isfolder(folder), mkdir(folder); end
file=fullfile(folder,sprintf('candidate_%04d_trial%02d.mat', ...
    double(row.candidate_id),double(bankrow.trial)));
if isfile(file), z=load(file,'est'); est=z.est; return; end
z=load(char(bankrow.file),'record');
est=fit_fixed_candidate(method,z.record,truth,row,cfg);
if cfg.save_fit_cache, save(file,'est','-v7.3'); end
end

function tau=tau_from_row(r)
tau=struct('supra',double(r.tau_supra),'within',double(r.tau_within), ...
    'cross',double(r.tau_cross),'copy',double(r.tau_copy));
end

function T=base_row(cfg,method,S,L,trial,tau)
T=table(cfg.anchor_p,method,double(S.candidate_id),string(S.selection_criterion), ...
    double(S.selection_rank),double(L.level),double(L.kappa_index),double(L.kappa), ...
    double(L.snr_index),double(L.target_snr_db),double(L.noise_amplitude_ratio), ...
    double(L.noise_to_signal_power_ratio),double(L.noise_energy_fraction), ...
    "viewwise_kappa_heteroscedastic_gaussian", ...
    double(L.unrotated_view_variance_1),double(L.unrotated_view_variance_2), ...
    double(L.unrotated_view_variance_3),double(L.unrotated_view_variance_4), ...
    "p5000_validation_fixed",double(trial),tau.supra,tau.within,tau.cross,tau.copy, ...
    'VariableNames',{'anchor_p','method','candidate_id','selection_criterion', ...
    'selection_rank','level','kappa_index','kappa','snr_index','target_snr_db', ...
    'noise_amplitude_ratio','noise_to_signal_power_ratio','noise_energy_fraction', ...
    'noise_model','unrotated_view_variance_1','unrotated_view_variance_2', ...
    'unrotated_view_variance_3','unrotated_view_variance_4','threshold_policy', ...
    'trial','tau_supra','tau_within','tau_cross','tau_copy'});
end

function T=tau_row(method,S,tau)
T=table(method,double(S.candidate_id),string(S.selection_criterion), ...
    double(S.selection_rank),"p5000_validation_fixed",tau.supra,tau.within, ...
    tau.cross,tau.copy,'VariableNames',{'method','candidate_id', ...
    'selection_criterion','selection_rank','threshold_policy','tau_supra', ...
    'tau_within','tau_cross','tau_copy'});
end

function T=add_realization_metrics(T,x)
n={'auc_supra','auc_within','auc_cross','auc_copy','auc_component', ...
   'f1_supra','f1_within','f1_cross','f1_copy','f1_component'};
for k=1:numel(n), T.(n{k})=x(k); end
end

function T=add_summary_metrics(T,X)
n={'auc_supra','auc_within','auc_cross','auc_copy','auc_component', ...
   'f1_supra','f1_within','f1_cross','f1_copy','f1_component'};
for k=1:numel(n)
    T.([n{k},'_mean'])=mean(X(:,k),'omitnan');
    T.([n{k},'_std'])=std(X(:,k),0,'omitnan');
end
end

function write_protocol(cfg,levels)
fid=fopen(fullfile(cfg.output_dir,'EXPERIMENT_PROTOCOL.txt'),'w');
if fid<0, return; end
c=onCleanup(@()fclose(fid));
fprintf(fid,'View-wise kappa x SNR heteroscedastic Gaussian-noise experiment\n');
fprintf(fid,'Signal count: p=%d; test realizations per condition: %d.\n',cfg.anchor_p,cfg.test_count);
fprintf(fid,'SNR grid (dB): %s\n',mat2str(cfg.snr_db));
fprintf(fid,'Kappa grid: %s\n',mat2str(cfg.kappa_grid));
fprintf(fid,'Endpoint view-variance profile v*: %s.\n',mat2str(cfg.view_variance_endpoint));
fprintf(fid,'Variance law: v(kappa)=1+kappa*(v*-1).\n');
fprintf(fid,'kappa=0 gives [1 1 1 1] (equal-variance IID/AWGN endpoint).\n');
fprintf(fid,'kappa=1 gives the full endpoint view-heteroscedastic profile.\n');
fprintf(fid,'Noise is VIEW-dependent only: all n nodes in one view share the same prescribed variance multiplier.\n');
fprintf(fid,'No node-wise low/high-noise partition is used.\n');
fprintf(fid,'Y_observed=Y_clean+alpha*E_kappa, alpha=10^(-SNR_dB/20).\n');
fprintf(fid,'The endpoint/interpolated view profile is cyclically rotated across trials: %d.\n',cfg.rotate_variance_profile);
fprintf(fid,'Within a trial, Y_clean and the underlying standard-Gaussian innovation E0 are paired across every kappa and SNR.\n');
fprintf(fid,'For each kappa, E_kappa is realization-wise normalized to the empirical power of Y_clean before SNR scaling.\n');
fprintf(fid,'Physical graph, B, Gamma, signal model, selected candidates, solver hyperparameters, and tau thresholds are frozen.\n');
fprintf(fid,'Configurations and thresholds are exactly those selected at p=5000 by mean validation component F1.\n');
fprintf(fid,'No hyperparameter selection or threshold recalibration occurs in this sweep.\n');
fprintf(fid,'Number of kappa x SNR conditions: %d.\n',height(levels));
clear c
end

function label=snr_label(value)
if isinf(value), label='clean'; else, label=sprintf('%g',value); end
end

function T=append_rows(T,S)
if isempty(T), T=S; else, T=[T;S]; end
end
