function summary=run_model_selection(cfg)
%RUN_MODEL_SELECTION Validation-based model selection at each sample size.
setup_three_method_paths();
if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end
[truth,index]=prepare_mosaic_fixed_truth_bank(cfg);
C=struct('MOSAIC',mosaic_planted_hyperparameter_catalog(cfg.mode), ...
    'PGL2021',pgl2021_broad_hyperparameter_catalog(cfg.mode), ...
    'ZW2024',zw2024_broad_hyperparameter_catalog(cfg.mode));
writetable(C.MOSAIC,fullfile(cfg.output_dir,'MOSAIC_CATALOG_USED.csv'));
writetable(C.PGL2021,fullfile(cfg.output_dir,'PGL2021_CATALOG_USED.csv'));
writetable(C.ZW2024,fullfile(cfg.output_dir,'ZW2024_CATALOG_USED.csv'));
save(fullfile(cfg.output_dir,'RUN_CONFIG_AND_TRUTH.mat'),'cfg','truth','-v7.3');
write_protocol(cfg,truth,C);
diary(fullfile(cfg.output_dir,'FULL_CONSOLE_LOG.txt'));
clean=onCleanup(@()diary('off')); %#ok<NASGU>
methods=["MOSAIC","PGL2021","ZW2024"]; alltest=table(); allsel=table();
for pi=1:numel(cfg.sample_counts)
    p=cfg.sample_counts(pi); pdir=fullfile(cfg.output_dir,sprintf('p_%05d',p));
    if ~exist(pdir,'dir'), mkdir(pdir); end
    val=index(index.sample_count==p&index.split=="validation",:);
    tst=index(index.sample_count==p&index.split=="test",:);
    screenbank=val(1:min(cfg.screen_validation_count,height(val)),:);
    fprintf('\n========== p=%d | validation=%d test=%d ==========\n',p,height(val),height(tst));
    for mi=1:numel(methods)
        method=methods(mi); catalog=C.(char(method));
        fprintf('\n--- %s screen: %d candidates ---\n',method,height(catalog));
        screen=evaluate_catalog(method,catalog,screenbank,truth,cfg,'screen',p,pdir);
        writetable(screen,fullfile(pdir,method+"_SCREEN_ALL.csv"));
        shortlist=shortlist_union(screen,cfg.shortlist_per_criterion);
        writetable(shortlist,fullfile(pdir,method+"_REFINEMENT_SHORTLIST.csv"));
        fprintf('--- %s refinement: %d candidates ---\n',method,height(shortlist));
        refine=evaluate_catalog(method,shortlist,val,truth,cfg,'refine',p,pdir);
        writetable(refine,fullfile(pdir,method+"_REFINEMENT_ALL.csv"));
        selected=select_component_f1_top3(refine,cfg.top_k);
        writetable(selected,fullfile(pdir,method+"_F1_COMPONENT_TOP3.csv"));
        allsel=append_rows(allsel,generic_selection_rows(selected,method,p));
        [TS,TR]=evaluate_selected_test(method,selected,tst,truth,cfg,p,pdir);
        writetable(TS,fullfile(pdir,method+"_TOP3_TEST_SUMMARY.csv"));
        writetable(TR,fullfile(pdir,method+"_TOP3_TEST_REALIZATIONS.csv"));
        alltest=append_rows(alltest,TS);
        if cfg.make_figures
            save_selected_recovery_figures(method,selected,tst,truth,cfg,p,pdir);
        end
    end
end
writetable(allsel,fullfile(cfg.output_dir,'ALL_F1_COMPONENT_TOP3_VALIDATION_SELECTIONS.csv'));
writetable(alltest,fullfile(cfg.output_dir,'ALL_F1_COMPONENT_TOP3_TEST_SUMMARY.csv'));
if cfg.make_figures, plot_three_method_score_curves(alltest,cfg.output_dir); end
cvx_report=zw_cvx_spotcheck(allsel,index,truth,cfg);
save(fullfile(cfg.output_dir,'COMPLETE_BENCHMARK.mat'), ...
    'cfg','truth','C','allsel','alltest','cvx_report','-v7.3');
fid=fopen(fullfile(cfg.output_dir,'BENCHMARK_COMPLETED.txt'),'w');
if fid>=0
    fprintf(fid,'Completed %s\n',char(datetime('now')));
    fprintf(fid,'Fixed realized topology: yes\nValidation/test leakage: no\n');
    fprintf(fid,'p values: %s\n',mat2str(cfg.sample_counts)); fclose(fid);
end
summary=struct('output_dir',cfg.output_dir,'test_summary',alltest, ...
    'validation_selections',allsel,'cvx_spotcheck',cvx_report);
end

function T=evaluate_catalog(method,catalog,bank,truth,cfg,stage,p,pdir)
file=fullfile(pdir,sprintf('%s_%s_CHECKPOINT.csv',method,upper(stage)));
if isfile(file), T=readtable(file,'TextType','string'); else, T=table(); end
done=[]; if ~isempty(T), done=double(T.candidate_id); end
for c=1:height(catalog)
    cid=double(catalog.candidate_id(c)); if ismember(cid,done), continue; end
    estimates=cell(height(bank),1); ok=false(height(bank),1);
    runtime=nan(height(bank),1); statuses=strings(height(bank),1);
    auc=nan(height(bank),5);
    for j=1:height(bank)
        try
            estimates{j}=load_or_fit(method,catalog(c,:),bank(j,:),truth,cfg,stage,p,pdir);
            ok(j)=true; runtime(j)=estimates{j}.runtime;
            statuses(j)=string(estimates{j}.status);
            a=graph_auc_metrics(estimates{j}.A_primary,truth);
            auc(j,:)=[a.supra,a.within,a.cross,a.copy,a.component_average];
        catch ME
            statuses(j)="error:"+string(ME.identifier);
            if cfg.verbose
                fprintf('  %s id=%d trial=%d failed: %s\n',method,cid,j,ME.message);
            end
        end
    end
    row=catalog(c,:); row.success_count=nnz(ok); row.validation_count=height(bank);
    row.eligible=all(ok);
    row.converged_fraction=mean(statuses=="converged"|statuses=="completed");
    row.mean_runtime_seconds=mean(runtime,'omitnan'); row.statuses=strjoin(statuses,'|');
    if all(ok)
        f=tune_graph_f1_thresholds(estimates,truth,cfg.threshold_grid);
        row.tau_supra=f.thresholds.supra; row.tau_within=f.thresholds.within;
        row.tau_cross=f.thresholds.cross; row.tau_copy=f.thresholds.copy;
        row.validation_f1_supra_mean=f.supra_mean;
        row.validation_f1_supra_std=f.supra_std;
        row.validation_f1_component_mean=f.component_average_mean;
        row.validation_f1_component_std=f.component_average_std;
        row.validation_auc_supra_mean=mean(auc(:,1),'omitnan');
        row.validation_auc_supra_std=std(auc(:,1),0,'omitnan');
        row.validation_auc_component_mean=mean(auc(:,5),'omitnan');
        row.validation_auc_component_std=std(auc(:,5),0,'omitnan');
        row.validation_auc_within_mean=mean(auc(:,2),'omitnan');
        row.validation_auc_cross_mean=mean(auc(:,3),'omitnan');
        row.validation_auc_copy_mean=mean(auc(:,4),'omitnan');
        row.validation_f1_within_mean=f.within_mean;
        row.validation_f1_cross_mean=f.cross_mean;
        row.validation_f1_copy_mean=f.copy_mean;
    else
        names={'tau_supra','tau_within','tau_cross','tau_copy', ...
            'validation_f1_supra_mean','validation_f1_supra_std', ...
            'validation_f1_component_mean','validation_f1_component_std', ...
            'validation_auc_supra_mean','validation_auc_supra_std', ...
            'validation_auc_component_mean','validation_auc_component_std', ...
            'validation_auc_within_mean','validation_auc_cross_mean', ...
            'validation_auc_copy_mean','validation_f1_within_mean', ...
            'validation_f1_cross_mean','validation_f1_copy_mean'};
        for k=1:numel(names), row.(names{k})=NaN; end
    end
    T=append_rows(T,row); writetable(T,file);
    if cfg.verbose
        fprintf('  %s %s id=%d | AUC %.3f/%.3f F1 %.3f/%.3f | %d/%d\n', ...
            method,stage,cid,row.validation_auc_supra_mean, ...
            row.validation_auc_component_mean,row.validation_f1_supra_mean, ...
            row.validation_f1_component_mean,row.success_count,row.validation_count);
    end
end
T=sortrows(T,'candidate_id');
end

function est=load_or_fit(method,row,bankrow,truth,cfg,stage,~,pdir)
folder=fullfile(pdir,'fit_cache',char(method),stage);
if ~exist(folder,'dir'), mkdir(folder); end
file=fullfile(folder,sprintf('candidate_%04d_trial_%02d.mat', ...
    double(row.candidate_id),double(bankrow.trial)));
if isfile(file), z=load(file,'est'); est=z.est; return; end
z=load(char(bankrow.file),'record');
est=fit_three_method_candidate(method,z.record,truth,row,cfg,stage);
save(file,'est','-v7.3');
end

function S=shortlist_union(T,k)
criteria={'validation_auc_supra_mean','validation_auc_component_mean', ...
    'validation_f1_supra_mean','validation_f1_component_mean'};
ids=[];
for j=1:numel(criteria)
    E=T(logical(T.eligible)&isfinite(T.(criteria{j})),:);
    if isempty(E), error('No eligible candidates for %s.',criteria{j}); end
    E=sortrows(E,{criteria{j},'converged_fraction','candidate_id'}, ...
        {'descend','descend','ascend'});
    ids=[ids;double(E.candidate_id(1:min(k,height(E))))]; %#ok<AGROW>
end
ids=unique(ids,'stable'); S=T(ismember(double(T.candidate_id),ids),:);
end

function S=select_component_f1_top3(T,k)
E=T(logical(T.eligible)&isfinite(T.validation_f1_component_mean),:);
if isempty(E), error('No eligible refined candidates for component F1.'); end
E=sortrows(E,{'validation_f1_component_mean','validation_f1_component_std', ...
    'converged_fraction','candidate_id'}, ...
    {'descend','ascend','descend','ascend'});
S=E(1:min(k,height(E)),:);
S.selection_criterion=repmat("f1_component",height(S),1);
S.selection_rank=(1:height(S))';
S.validation_selection_score=S.validation_f1_component_mean;
end

function [TS,TR]=evaluate_selected_test(method,S,bank,truth,cfg,p,pdir)
TS=table(); TR=table();
for q=1:height(S)
    tau=struct('supra',S.tau_supra(q),'within',S.tau_within(q), ...
        'cross',S.tau_cross(q),'copy',S.tau_copy(q));
    scores=nan(height(bank),10); runtimes=nan(height(bank),1);
    statuses=strings(height(bank),1);
    for j=1:height(bank)
        try
            est=load_or_fit(method,S(q,:),bank(j,:),truth,cfg,'test',p,pdir);
            a=graph_auc_metrics(est.A_primary,truth);
            f=apply_graph_f1_thresholds(est.A_primary,truth,tau);
            scores(j,:)=[a.supra,a.within,a.cross,a.copy,a.component_average, ...
                f.supra,f.within,f.cross,f.copy,f.component_average];
            runtimes(j)=est.runtime; statuses(j)=string(est.status);
        catch ME
            statuses(j)="error:"+string(ME.identifier);
        end
        rr=table(method,p,double(S.candidate_id(q)),S.selection_criterion(q), ...
            double(S.selection_rank(q)),double(bank.trial(j)), ...
            scores(j,1),scores(j,2),scores(j,3),scores(j,4),scores(j,5), ...
            scores(j,6),scores(j,7),scores(j,8),scores(j,9),scores(j,10), ...
            'VariableNames',{'method','sample_count','candidate_id','selection_criterion', ...
            'selection_rank','trial','auc_supra','auc_within','auc_cross','auc_copy', ...
            'auc_component','f1_supra','f1_within','f1_cross','f1_copy','f1_component'});
        TR=append_rows(TR,rr);
    end
    vals={method,p,double(S.candidate_id(q)),S.selection_criterion(q), ...
        double(S.selection_rank(q)),double(S.validation_selection_score(q)), ...
        tau.supra,tau.within,tau.cross,tau.copy};
    names={'method','sample_count','candidate_id','selection_criterion','selection_rank', ...
        'validation_selection_score','tau_supra','tau_within','tau_cross','tau_copy'};
    metric={'auc_supra','auc_within','auc_cross','auc_copy','auc_component', ...
        'f1_supra','f1_within','f1_cross','f1_copy','f1_component'};
    for k=1:10
        vals=[vals,{mean(scores(:,k),'omitnan'),std(scores(:,k),0,'omitnan')}]; %#ok<AGROW>
        names=[names,{['test_',metric{k},'_mean'],['test_',metric{k},'_std']}]; %#ok<AGROW>
    end
    vals=[vals,{mean(runtimes,'omitnan'),strjoin(statuses,'|')}];
    names=[names,{'mean_runtime_seconds','statuses'}];
    TS=append_rows(TS,cell2table(vals,'VariableNames',names));
end
end

function G=generic_selection_rows(S,method,p)
G=table(repmat(method,height(S),1),repmat(p,height(S),1),double(S.candidate_id), ...
    S.selection_criterion,double(S.selection_rank),double(S.validation_selection_score), ...
    S.tau_supra,S.tau_within,S.tau_cross,S.tau_copy, ...
    'VariableNames',{'method','sample_count','candidate_id','selection_criterion', ...
    'selection_rank','validation_selection_score','tau_supra','tau_within','tau_cross','tau_copy'});
end

function report=zw_cvx_spotcheck(S,index,truth,cfg)
report=struct('attempted',false,'available',exist('cvx_begin','file')==2);
if ~cfg.cvx_spotcheck||~report.available, return; end
Z=S(S.method=="ZW2024"&S.sample_count==max(cfg.sample_counts)& ...
    S.selection_criterion=="f1_component"&S.selection_rank==1,:);
if isempty(Z), return; end
cat=zw2024_broad_hyperparameter_catalog(cfg.mode);
row=cat(cat.candidate_id==Z.candidate_id(1),:);
B=index(index.split=="validation"&index.sample_count==max(cfg.sample_counts),:);
z=load(char(B.file(1)),'record');
fast=fit_three_method_candidate("ZW2024",z.record,truth,row,cfg,'refine');
c2=cfg; c2.zw_backend='cvx'; cvxfit=fit_three_method_candidate("ZW2024",z.record,truth,row,c2,'refine');
report.attempted=true;
report.relative_frobenius=norm(fast.A_primary-cvxfit.A_primary,'fro')/ ...
    max(1,norm(cvxfit.A_primary,'fro'));
end

function write_protocol(cfg,t,C)
fid=fopen(fullfile(cfg.output_dir,'EXPERIMENT_PROTOCOL.txt'),'w'); if fid<0, return; end
fprintf(fid,'MOSAIC realized graph: %s\n',cfg.source_dataset);
fprintf(fid,'n=%d K=%d r=%d N=%d\n',t.n,t.K,t.r,t.N);
fprintf(fid,'p=%s; validation/test=%d/%d; nested prefixes within realization.\n', ...
    mat2str(cfg.sample_counts),cfg.validation_count,cfg.test_count);
fprintf(fid,'Catalog sizes MOSAIC/PGL2021/ZW2024: %d/%d/%d\n', ...
    height(C.MOSAIC),height(C.PGL2021),height(C.ZW2024));
fprintf(fid,'The top three candidates minimize the validation selection risk through mean component F1.\n');
fprintf(fid,'Component-specific relative thresholds are selected on validation data and applied to test data.\n');
fprintf(fid,'ZW2024 primary score is its effective interaction matrix.\n'); fclose(fid);
end
function T=append_rows(T,S), if isempty(T), T=S; else, T=[T;S]; end, end
