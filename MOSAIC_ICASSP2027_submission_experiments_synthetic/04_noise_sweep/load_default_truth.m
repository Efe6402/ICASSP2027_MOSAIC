function truth=load_default_truth(cfg)
z=load(cfg.default_truth_file);
if isfield(z,'truth'), d=z.truth; elseif isfield(z,'data'), d=z.data;
else, error('Default truth file must contain truth or data.'); end
need={'n','K','r','A_state_true','B_true','Gamma_true','A_supra_true', ...
    'A_physical_true','A_copy_true','truth_masks','truth_support_structural', ...
    'Sigma','Sigma_clean','cfg','community_nodes'};
for j=1:numel(need), assert(isfield(d,need{j}),'Missing truth field %s.',need{j}); end
truth=d; truth.N=truth.n*truth.K;
truth.A_within_true=truth.A_supra_true.*truth.truth_masks.within;
truth.A_cross_true=truth.A_supra_true.*truth.truth_masks.crossnode_crossview;
truth.A_copy_masked_true=truth.A_supra_true.*truth.truth_masks.copy;
end
