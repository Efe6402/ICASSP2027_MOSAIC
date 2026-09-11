function truth=build_rho_truth(base,rho)
%BUILD_RHO_TRUTH Rebuild B, supra-adjacency, and covariance for one rho.
assert(base.K==4&&base.r==2, ...
    'This controlled profile path is defined for K=4 and r=2.');
assert(rho>=0&&rho<=.5,'rho must lie in [0,0.5].');
W=[1,0;1-rho,rho;rho,1-rho;0,1];
B=sqrt(W); B=B./sum(B,1);
is_reference=abs(rho-5/102)<1e-12;
if is_reference, B=base.B_true; end

truth=base; truth.rho=double(rho); truth.W_true=W; truth.B_true=B;
truth.generator_tag=sprintf('B_profile_rho_%.12g',rho);
n=truth.n; K=truth.K; r=truth.r; N=n*K;
Aphysical=zeros(N); Acopy=zeros(N);
for m=1:r
    Q=B(:,m)*B(:,m)';
    Aphysical=Aphysical+kron(Q,truth.A_state_true(:,:,m));
    Acopy=Acopy+kron(Q-diag(diag(Q)),diag(truth.Gamma_true(:,m)));
end
truth.A_physical_true=zero_sym(Aphysical);
truth.A_copy_true=zero_sym(Acopy);
truth.A_supra_true=zero_sym(Aphysical+Acopy);
if is_reference
    % The reference level uses the embedded matrices directly.
    truth.A_physical_true=base.A_physical_true;
    truth.A_copy_true=base.A_copy_true;
    truth.A_supra_true=base.A_supra_true;
end
truth.truth_support_structural=struct( ...
    'all',truth.A_supra_true>0, ...
    'within',(truth.A_supra_true>0)&truth.truth_masks.within, ...
    'crossnode_crossview',(truth.A_supra_true>0)&truth.truth_masks.crossnode_crossview, ...
    'copy',(truth.A_supra_true>0)&truth.truth_masks.copy);
truth.A_within_true=truth.A_supra_true.*truth.truth_masks.within;
truth.A_cross_true=truth.A_supra_true.*truth.truth_masks.crossnode_crossview;
truth.A_copy_masked_true=truth.A_supra_true.*truth.truth_masks.copy;
truth.N=N;
if is_reference
    truth.Sigma_clean=base.Sigma_clean; truth.Sigma=base.Sigma;
else
    [truth.Sigma_clean,truth.Sigma]=signal_covariance_from_supra( ...
        truth.A_supra_true,truth.cfg);
end
end

function A=zero_sym(A)
A=.5*(A+A'); A(1:size(A,1)+1:end)=0;
end
