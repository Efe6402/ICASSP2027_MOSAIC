function m=graph_auc_metrics(A,truth)
%GRAPH_AUC_METRICS Tie-aware AUC on strict-upper component masks.
A=max(0,.5*(A+A')); N=size(A,1); upper=triu(true(N),1);
defs={'supra','all','all';'within','within','within'; ...
    'cross','crossnode_crossview','crossnode_crossview';'copy','copy','copy'};
m=struct();
for j=1:size(defs,1)
    mask=upper&truth.truth_masks.(defs{j,3});
    y=logical(truth.truth_support_structural.(defs{j,2})(mask));
    m.(defs{j,1})=rank_auc(A(mask),y);
end
m.component=mean([m.within,m.cross,m.copy],'omitnan');
end
function a=rank_auc(x,y)
x=x(:); y=logical(y(:)); np=nnz(y); nn=nnz(~y);
if np==0||nn==0, a=NaN; return; end
[xs,ord]=sort(x,'ascend'); ranks=zeros(size(x)); i=1;
while i<=numel(x)
    j=i; while j<numel(x)&&xs(j+1)==xs(i), j=j+1; end
    ranks(ord(i:j))=(i+j)/2; i=j+1;
end
a=(sum(ranks(y))-np*(np+1)/2)/(np*nn);
end
