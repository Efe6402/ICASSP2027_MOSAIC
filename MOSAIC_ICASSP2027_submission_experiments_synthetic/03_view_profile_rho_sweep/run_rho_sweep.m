function summary=run_rho_sweep(cfg)
%RUN_RHO_SWEEP Evaluate p=5000-selected models along the rho path.
setup_rho_sweep_paths(cfg.root);
if ~isfolder(cfg.output_dir), mkdir(cfg.output_dir); end
[truths,index,levels]=prepare_rho_signal_bank(cfg);
bank_validation=validate_rho_bank(cfg,truths,index,levels);
writetable(levels,fullfile(cfg.output_dir,'RHO_LEVELS.csv'));
writetable(index,fullfile(cfg.output_dir,'BANK_INDEX.csv'));
writetable(struct2table(bank_validation),fullfile(cfg.output_dir,'BANK_VALIDATION.csv'));
save(fullfile(cfg.output_dir,'RUN_CONFIG.mat'),'cfg','-v7.3');
if cfg.make_truth_figures, save_rho_truth_figures(cfg); end
diary(fullfile(cfg.output_dir,'FULL_CONSOLE_LOG.txt'));
cleaner=onCleanup(@()diary('off'));
methods=["MOSAIC","PGL2021","ZW2024"];
allsummary=table(); allreal=table(); tauaudit=table();
for method=methods
    selected=load_selected(cfg,method);
    writetable(selected,fullfile(cfg.output_dir,method+"_CONFIGURATIONS.csv"));
    for q=1:height(selected)
        tau=tau_from_row(selected(q,:));
        for lev=1:height(levels)
            truth=truths{lev}; tst=index(index.level==lev&index.split=="test",:);
            X=nan(height(tst),10); runtime=nan(height(tst),1); status=strings(height(tst),1);
            for j=1:height(tst)
                est=load_or_fit(cfg,method,selected(q,:),tst(j,:),truth,lev);
                a=graph_auc_metrics(est.A_primary,truth);
                f=apply_graph_f1_thresholds(est.A_primary,truth,tau);
                X(j,:)=[a.supra,a.within,a.cross,a.copy,a.component, ...
                    f.supra,f.within,f.cross,f.copy,f.component];
                runtime(j)=est.runtime; status(j)=string(est.status);
                row=base_row(cfg,method,selected(q,:),lev,levels.rho(lev),j,tau);
                row=add_realization_metrics(row,X(j,:));
                row.runtime_seconds=runtime(j); row.status=status(j);
                allreal=append_rows(allreal,row);
            end
            row=base_row(cfg,method,selected(q,:),lev,levels.rho(lev),NaN,tau);
            row=add_summary_metrics(row,X);
            row.mean_runtime_seconds=mean(runtime,'omitnan');
            row.converged_fraction=mean(status=="converged"|status=="completed");
            allsummary=append_rows(allsummary,row);
            tauaudit=append_rows(tauaudit,tau_row(cfg,method,selected(q,:),lev,levels.rho(lev),tau));
            writetable(allsummary,fullfile(cfg.output_dir,'TEST_SUMMARY_CHECKPOINT.csv'));
            writetable(allreal,fullfile(cfg.output_dir,'TEST_REALIZATIONS_CHECKPOINT.csv'));
        end
    end
end
writetable(allsummary,fullfile(cfg.output_dir,'TEST_SUMMARY.csv'));
writetable(allreal,fullfile(cfg.output_dir,'TEST_REALIZATIONS.csv'));
writetable(tauaudit,fullfile(cfg.output_dir,'TAU_AUDIT.csv'));
if cfg.make_figures, plot_rho_results(allsummary,cfg,cfg.output_dir); end
save(fullfile(cfg.output_dir,'COMPLETE_RESULTS.mat'),'cfg','levels','bank_validation', ...
    'allsummary','allreal','tauaudit','-v7.3');
write_protocol(cfg);
summary=struct('output_dir',cfg.output_dir,'rho_levels',levels, ...
    'test_summary',allsummary,'test_realizations',allreal,'thresholds',tauaudit);
end

function C=load_selected(cfg,method)
file=fullfile(cfg.selected_config_dir,method+"_F1_COMPONENT_TOP3.csv");
assert(isfile(file),'Missing selected-configuration table: %s',file);
C=readtable(file,'TextType','string');
C=C(C.selection_criterion=="f1_component"&C.selection_rank<=cfg.top_k,:);
C=sortrows(C,'selection_rank'); assert(height(C)==cfg.top_k);
end

function est=load_or_fit(cfg,method,row,bankrow,truth,lev)
folder=fullfile(cfg.output_dir,'FIT_CACHE',char(method),sprintf('rho_%02d',lev));
if ~isfolder(folder), mkdir(folder); end
file=fullfile(folder,sprintf('candidate_%04d_trial%02d.mat', ...
    double(row.candidate_id),double(bankrow.trial)));
if isfile(file), z=load(file,'est'); est=z.est; return; end
z=load(char(bankrow.file),'record'); est=fit_fixed_candidate(method,z.record,truth,row,cfg);
if cfg.save_fit_cache, save(file,'est','-v7.3'); end
end

function tau=tau_from_row(r)
tau=struct('supra',double(r.tau_supra),'within',double(r.tau_within), ...
    'cross',double(r.tau_cross),'copy',double(r.tau_copy));
end

function T=base_row(cfg,method,S,lev,rho,trial,tau)
T=table(cfg.anchor_p,method,double(S.candidate_id),string(S.selection_criterion), ...
    double(S.selection_rank),lev,double(rho),"p5000_validation_fixed",double(trial), ...
    tau.supra,tau.within,tau.cross,tau.copy,'VariableNames', ...
    {'anchor_p','method','candidate_id','selection_criterion','selection_rank', ...
    'level','rho','threshold_policy','trial','tau_supra','tau_within','tau_cross','tau_copy'});
end

function T=tau_row(cfg,method,S,lev,rho,tau)
T=base_row(cfg,method,S,lev,rho,NaN,tau);
T(:,{'trial'})=[];
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

function write_protocol(cfg)
fid=fopen(fullfile(cfg.output_dir,'EXPERIMENT_PROTOCOL.txt'),'w'); if fid<0, return; end
c=onCleanup(@()fclose(fid));
fprintf(fid,'Signal count p=%d; validation/test=%d/%d.\n', ...
    cfg.anchor_p,cfg.validation_count,cfg.test_count);
fprintf(fid,'rho values: %s\n',mat2str(cfg.rho_values));
fprintf(fid,'Physical modal graphs and copy coefficients are held fixed.\n');
fprintf(fid,'Ground-truth view profiles and the resulting supra-adjacency vary with rho.\n');
fprintf(fid,'Candidate identities, hyperparameters, and thresholds are selected at p=5000.\n');
fprintf(fid,'Selection criterion: mean validation component F1.\n');
end

function T=append_rows(T,S), if isempty(T), T=S; else, T=[T;S]; end, end
