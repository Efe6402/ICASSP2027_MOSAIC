function record=make_signal_record(Y,truth,meta)
p=size(Y,2); Cy=(Y*Y')/p; d=diag(Cy); S=max(d+d'-2*Cy,0);
S=.5*(S+S'); S(1:size(S,1)+1:end)=0;
Sp=zeros(truth.n); Sq=zeros(truth.K);
for s=1:p
    Xi=reshape(Y(:,s),truth.n,truth.K);
    Sp=Sp+cov(Xi'); Sq=Sq+cov(Xi);
end
record=meta; record.sample_count=p; record.n=truth.n; record.K=truth.K;
record.S=S; record.Cy_uncentered=Cy; record.Sp=Sp/p; record.Sq=Sq/p;
end
