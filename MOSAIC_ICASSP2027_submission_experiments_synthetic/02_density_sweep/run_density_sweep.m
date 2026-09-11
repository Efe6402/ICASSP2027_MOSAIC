function summary=run_density_sweep(cfg)
%RUN_DENSITY_SWEEP Evaluate p=5000-selected models along a density path.
setup_three_method_paths();
if ~isfolder(cfg.output_dir), mkdir(cfg.output_dir); end
[truths,index,D]=prepare_density_bank(cfg);
writetable(D,fullfile(cfg.output_dir,'DENSITY_LEVELS.csv'));
save(fullfile(cfg.output_dir,'RUN_CONFIG.mat'),'cfg','-v7.3');
diary(fullfile(cfg.output_dir,'FULL_CONSOLE_LOG.txt'));
cleaner=onCleanup(@()diary('off'));
methods=["MOSAIC","PGL2021","ZW2024"];
allsummary=table(); allreal=table();
for d=1:numel(truths)
    truth=truths{d}; tst=index(index.density_id==d&index.split=="test",:);
    ddir=fullfile(cfg.output_dir,sprintf('density_%02d',d)); mkdir_if_needed(ddir);
    for method=methods
        selected=load_selected(cfg,method);
        writetable(selected,fullfile(ddir,method+"_CONFIGURATIONS.csv"));
        for q=1:height(selected)
            tau=struct('supra',double(selected.tau_supra(q)), ...
                'within',double(selected.tau_within(q)), ...
                'cross',double(selected.tau_cross(q)), ...
                'copy',double(selected.tau_copy(q)));
            X=nan(height(tst),10); runtime=nan(height(tst),1); status=strings(height(tst),1);
            for j=1:height(tst)
                est=load_or_fit(method,selected(q,:),tst(j,:),truth,cfg,d,ddir);
                a=graph_auc_metrics(est.A_primary,truth);
                f=apply_graph_f1_thresholds(est.A_primary,truth,tau);
                X(j,:)=[a.supra,a.within,a.cross,a.copy,a.component_average, ...
                    f.supra,f.within,f.cross,f.copy,f.component_average];
                runtime(j)=est.runtime; status(j)=string(est.status);
                row=base_row(method,selected(q,:),D(d,:),j,tau);
                row=add_realization_metrics(row,X(j,:));
                row.runtime_seconds=runtime(j); row.status=status(j);
                allreal=append_rows(allreal,row);
            end
            row=base_row(method,selected(q,:),D(d,:),NaN,tau);
            row=add_summary_metrics(row,X);
            row.mean_runtime_seconds=mean(runtime,'omitnan');
            row.converged_fraction=mean(status=="converged"|status=="completed");
            allsummary=append_rows(allsummary,row);
            writetable(allsummary,fullfile(cfg.output_dir,'TEST_SUMMARY_CHECKPOINT.csv'));
            writetable(allreal,fullfile(cfg.output_dir,'TEST_REALIZATIONS_CHECKPOINT.csv'));
        end
    end
end
writetable(allsummary,fullfile(cfg.output_dir,'TEST_SUMMARY.csv'));
writetable(allreal,fullfile(cfg.output_dir,'TEST_REALIZATIONS.csv'));
if cfg.make_figures, plot_density_results(allsummary,cfg.output_dir); end
save(fullfile(cfg.output_dir,'COMPLETE_RESULTS.mat'),'cfg','D','allsummary','allreal','-v7.3');
write_protocol(cfg,D);
summary=struct('output_dir',cfg.output_dir,'density_levels',D, ...
    'test_summary',allsummary,'test_realizations',allreal);
end

function C=load_selected(cfg,method)
file=fullfile(cfg.selected_config_dir,method+"_F1_COMPONENT_TOP3.csv");
assert(isfile(file),'Missing selected-configuration table: %s',file);
C=readtable(file,'TextType','string');
C=C(C.selection_criterion=="f1_component"&C.selection_rank<=cfg.top_k,:);
C=sortrows(C,'selection_rank'); assert(height(C)==cfg.top_k);
end

function est=load_or_fit(method,row,bankrow,truth,cfg,d,ddir)
folder=fullfile(ddir,'fit_cache',char(method)); mkdir_if_needed(folder);
file=fullfile(folder,sprintf('candidate_%04d_trial_%02d.mat', ...
    double(row.candidate_id),double(bankrow.trial)));
if isfile(file), z=load(file,'est'); est=z.est; return; end
z=load(char(bankrow.file),'record');
est=fit_three_method_candidate(method,z.record,truth,row,cfg,'fixed_density');
est.density_id=d; save(file,'est','-v7.3');
end

function T=base_row(method,S,D,trial,tau)
T=table("fixed", "p5000_validation_fixed",method,double(D.density_id), ...
    double(D.q_in_mode1),double(D.q_in_mode2),double(D.q_out), ...
    double(D.expected_mean_modal_density),double(D.realized_mean_modal_density), ...
    double(S.candidate_id),string(S.selection_criterion),double(S.selection_rank), ...
    double(trial),tau.supra,tau.within,tau.cross,tau.copy, ...
    'VariableNames',{'phase','threshold_policy','method','density_id', ...
    'q_in_mode1','q_in_mode2','q_out','expected_density','realized_density', ...
    'candidate_id','selection_criterion','selection_rank','trial', ...
    'tau_supra','tau_within','tau_cross','tau_copy'});
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
    T.(['test_',n{k},'_mean'])=mean(X(:,k),'omitnan');
    T.(['test_',n{k},'_std'])=std(X(:,k),0,'omitnan');
end
end

function write_protocol(cfg,D)
fid=fopen(fullfile(cfg.output_dir,'EXPERIMENT_PROTOCOL.txt'),'w'); if fid<0, return; end
c=onCleanup(@()fclose(fid));
fprintf(fid,'Signal count p=%d; test realizations=%d.\n',cfg.signal_count,cfg.test_count);
fprintf(fid,'Modal edge-count schedule: %s\n',mat2str(cfg.target_mode_edges));
fprintf(fid,'Block/background odds ratio: %.8g.\n',cfg.contrast_odds_ratio);
fprintf(fid,'Ground-truth copy coefficients and view profiles are held fixed.\n');
fprintf(fid,'Candidate identities, hyperparameters, and thresholds are selected at p=5000.\n');
fprintf(fid,'Selection criterion: mean validation component F1.\n');
fprintf(fid,'Realized densities: %s\n',mat2str(D.realized_mean_modal_density',5));
end

function mkdir_if_needed(p), if ~isfolder(p), mkdir(p); end, end
function T=append_rows(T,S), if isempty(T), T=S; else, T=[T;S]; end, end
