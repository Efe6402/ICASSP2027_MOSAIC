function T=pgl2021_broad_hyperparameter_catalog(mode)
%PGL2021_BROAD_HYPERPARAMETER_CATALOG Original region plus off-domain grid.
if nargin<1, mode='full'; end
wide=logspace(-4,1.5,25);
fine1=.10:.02:.50; fine2=.10:.02:.40;
[A,B]=ndgrid(wide,wide); V=[A(:),B(:)];
[A,B]=ndgrid(fine1,fine2); V=[V;A(:),B(:);.20,.30;.25,.25];
V=unique(round(V,12),'rows','stable');
label=repmat("wide_log_and_authors_grid",size(V,1),1);
label(abs(V(:,1)-.2)<1e-12&abs(V(:,2)-.3)<1e-12)="paper_setting";
label(abs(V(:,1)-.25)<1e-12&abs(V(:,2)-.25)<1e-12)="repository_setting";
T=table((1:size(V,1))',label,V(:,1),V(:,2), ...
    'VariableNames',{'candidate_id','profile_name','beta1','beta2'});
if strcmpi(mode,'quick')
    keep=find(T.profile_name~="wide_log_and_authors_grid");
    keep=unique([keep;1;ceil(height(T)/2)]); T=T(keep,:);
    T.candidate_id=(1:height(T))';
end
end
