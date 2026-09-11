function T=zw2024_broad_hyperparameter_catalog(mode)
%ZW2024_BROAD_HYPERPARAMETER_CATALOG Broad AO search for off-domain data.
if nargin<1, mode='full'; end
lambda=[0 .01 .03 .1 .3 1 3 5 10 30];
rhoL=[1e-3 3e-3 1e-2 3e-2 .1 .3 1 3 10 30 100];
ratio=[1 3 5 10 15 30 100];
[L,R,Q]=ndgrid(lambda,rhoL,ratio);
V=[L(:),R(:),R(:).*Q(:)];
V=[V;.1,1,15;5,1,15;0,1,15];
V=unique(round(V,12),'rows','stable');
label=repmat("broad_lambda_rho_ratio",size(V,1),1);
label(abs(V(:,1)-.1)<1e-12&abs(V(:,2)-1)<1e-12&abs(V(:,3)-15)<1e-12)="paper_weak";
label(abs(V(:,1)-5)<1e-12&abs(V(:,2)-1)<1e-12&abs(V(:,3)-15)<1e-12)="paper_strong";
label(V(:,1)==0&abs(V(:,2)-1)<1e-12&abs(V(:,3)-15)<1e-12)="paper_lambda_zero";
T=table((1:size(V,1))',label,V(:,1),V(:,2),V(:,3), ...
    'VariableNames',{'candidate_id','profile_name','lambda','rho_L','rho_C'});
if strcmpi(mode,'quick')
    keep=T.profile_name~="broad_lambda_rho_ratio";
    T=T(keep,:); T.candidate_id=(1:height(T))';
end
end
