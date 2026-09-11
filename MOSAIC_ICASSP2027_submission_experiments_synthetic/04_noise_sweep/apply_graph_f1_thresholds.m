function m=apply_graph_f1_thresholds(A,truth,tau)
defs={'supra','all','all';'within','within','within'; ...
    'cross','crossnode_crossview','crossnode_crossview';'copy','copy','copy'};
m=struct();
for j=1:size(defs,1)
    z=one(A,truth,defs{j,2},defs{j,3},tau.(defs{j,1}));
    m.(defs{j,1})=z.f1; m.([defs{j,1},'_precision'])=z.precision;
    m.([defs{j,1},'_recall'])=z.recall;
end
m.component=mean([m.within,m.cross,m.copy]);
end
function m=one(A,t,sname,mname,tau)
N=size(A,1); mask=triu(true(N),1)&t.truth_masks.(mname);
x=max(0,A(mask)); mx=max(x); if mx>0, x=x/mx; else, x=zeros(size(x)); end
y=logical(t.truth_support_structural.(sname)(mask)); p=x>tau;
tp=nnz(p&y); fp=nnz(p&~y); fn=nnz(~p&y);
if tp==0, pr=0; re=0; else, pr=tp/(tp+fp); re=tp/(tp+fn); end
if pr+re==0, f=0; else, f=2*pr*re/(pr+re); end
m=struct('f1',f,'precision',pr,'recall',re);
end
