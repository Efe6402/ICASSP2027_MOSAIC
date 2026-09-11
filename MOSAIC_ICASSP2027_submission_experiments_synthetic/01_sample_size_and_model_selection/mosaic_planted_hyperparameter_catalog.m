function T=mosaic_planted_hyperparameter_catalog(mode)
%MOSAIC_PLANTED_HYPERPARAMETER_CATALOG Results-guided, truth-independent grid.
if nargin<1, mode='full'; end
names={'alpha_a','alpha_c','beta_a','beta_c','tau_intra','tau_inter', ...
    'lambda_B','eta_c','epsilon_intra','epsilon_inter','theta_scale','gamma_scale'};
w298=[.191302326773625,.00279057959231506,10.6600383849336,.0059175594578156, ...
    .0005654468572754,.0010795694734749,.318995274076651,1.67466537304686, ...
    .0010750399231159,.0012425138447346,.878512787093594,.0062952429896908];
w327=w298; w327(2)=9.30193197438352e-5; w327(12)=.062952429896908;
w296=w298; w296(2)=.000279057959231506;
balanced=[.209256854316246,.0073359352758823,.0524513366850577,.209037662441851, ...
    .0123483173496519,.270479880832943,.543692893240052,1.00649214234789, ...
    1.8039859684388e-7,2.772915731464e-4,.383771600299731,.302957527108981];
copy=[.01,.03,.1,.3,.1,.3,.01,1,1e-5,1e-5,.2,.2];
base=[w298;w327;w296;balanced;copy];
labels=["reference_candidate_298";"reference_candidate_327";"reference_candidate_296"; ...
    "balanced_anchor";"copy_anchor"];
lo=log10([1e-3,1e-5,1e-3,1e-4,1e-5,1e-5,1e-5,.1,1e-8,1e-8,.01,1e-3]);
hi=log10([1,.3,30,3,1,5,3,5,1e-2,1e-2,2,2]);
span=[.65,1.1,.65,1.1,1,1.3,.9,.35,.8,.9,.45,1.2];
V=base; source=labels;
for j=1:size(base,1)
    nlocal=48; if j<=3, nlocal=64; end
    X=local_design(nlocal,base(j,:),span,lo,hi,73000+j);
    V=[V;X]; source=[source;repmat("local_"+labels(j),nlocal,1)]; %#ok<AGROW>
end
% Geometric bridges between the supra winners and copy/balanced regimes.
for t=[.25 .5 .75]
    anchor=10.^((1-t)*log10(w298)+t*log10(copy));
    X=local_design(16,anchor,.45*ones(1,12),lo,hi,74000+round(100*t));
    V=[V;X]; source=[source;repmat("winner_copy_bridge",16,1)]; %#ok<AGROW>
end
X=stratified(64,lo,hi,75001); V=[V;X];
source=[source;repmat("global_safeguard",64,1)];
% Deterministic one-factor paths in the most sensitive copy/inter-view terms.
for idx=[2 4 6 7 8 10 12]
    for mult=[.1 .3 3 10 30]
        x=w298; x(idx)=min(max(x(idx)*mult,10^lo(idx)),10^hi(idx));
        V=[V;x]; source=[source;"one_factor_path"]; %#ok<AGROW>
    end
end
[~,keep]=unique(round(log10(V),10),'rows','stable'); V=V(keep,:); source=source(keep);
T=table((1:size(V,1))',source,'VariableNames',{'candidate_id','profile_name'});
for j=1:numel(names), T.(names{j})=V(:,j); end
T.r=repmat(2,height(T),1);
if strcmpi(mode,'quick')
    T=T(1:min(6,height(T)),:); T.candidate_id=(1:height(T))';
end
end
function X=local_design(N,b,span,lo,hi,seed)
X=stratified(N,max(log10(b)-span,lo),min(log10(b)+span,hi),seed);
end
function V=stratified(N,lo,hi,seed)
old=rng; c=onCleanup(@()rng(old)); rng(seed,'twister');
U=zeros(N,numel(lo));
for j=1:numel(lo), z=((0:N-1)'+rand(N,1))/N; U(:,j)=z(randperm(N)); end
V=10.^(lo+(hi-lo).*U);
end
