function out=tune_graph_f1_thresholds(estimates,truth,grid)
%TUNE_GRAPH_F1_THRESHOLDS Select one relative tau per component on validation.
defs={'supra','all','all';'within','within','within'; ...
    'cross','crossnode_crossview','crossnode_crossview';'copy','copy','copy'};
R=numel(estimates); out=struct('thresholds',struct(),'per_realization',struct());
for j=1:size(defs,1)
    label=defs{j,1}; scores=nan(R,numel(grid));
    for r=1:R
        for q=1:numel(grid)
            z=one_f1(estimates{r}.A_primary,truth,defs{j,2},defs{j,3},grid(q));
            scores(r,q)=z.f1;
        end
    end
    mu=mean(scores,1,'omitnan'); sd=std(scores,0,1,'omitnan');
    best=find(mu>=max(mu)-1e-14);
    [~,k]=min(sd(best)); best=best(k);
    % If mean and spread tie, prefer the more conservative larger threshold.
    tied=find(abs(mu-mu(best))<=1e-14 & abs(sd-sd(best))<=1e-14);
    best=tied(end);
    out.thresholds.(label)=grid(best);
    out.per_realization.(label)=scores(:,best);
    out.([label,'_mean'])=mu(best); out.([label,'_std'])=sd(best);
end
out.component_average_per=mean([out.per_realization.within, ...
    out.per_realization.cross,out.per_realization.copy],2);
out.component_average_mean=mean(out.component_average_per,'omitnan');
out.component_average_std=std(out.component_average_per,0,'omitnan');
end

function m=one_f1(A,truth,support_name,mask_name,tau)
N=size(A,1); upper=triu(true(N),1); mask=upper&truth.truth_masks.(mask_name);
x=max(0,A(mask)); mx=max(x); if mx>0, x=x/mx; else, x=zeros(size(x)); end
y=logical(truth.truth_support_structural.(support_name)(mask)); pred=x>tau;
tp=nnz(pred&y); fp=nnz(pred&~y); fn=nnz(~pred&y);
if tp==0, precision=0; recall=0; else, precision=tp/(tp+fp); recall=tp/(tp+fn); end
if precision+recall==0, f1=0; else, f1=2*precision*recall/(precision+recall); end
m=struct('f1',f1,'precision',precision,'recall',recall,'tp',tp,'fp',fp,'fn',fn);
end
