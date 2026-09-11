function report=validate_rho_bank(cfg,truths,index,levels)
%VALIDATE_RHO_BANK Fail fast on dimensions, seeds, profiles, and files.
assert(numel(truths)==numel(cfg.rho_values));
assert(height(levels)==numel(cfg.rho_values));
assert(height(index)==numel(cfg.rho_values)*(cfg.validation_count+cfg.test_count));
assert(all(index.sample_count==cfg.anchor_p));
assert(all(isfile(index.file)),'At least one bank record is missing.');
base=load_default_truth(cfg);
for lev=1:numel(truths)
    t=truths{lev}; assert(abs(t.rho-cfg.rho_values(lev))<1e-12);
    assert(all(abs(sum(t.B_true,1)-1)<1e-12),'B columns are not simplex normalized.');
    assert(isequal(t.A_state_true,base.A_state_true),'Physical modal graphs changed.');
    assert(isequal(t.Gamma_true,base.Gamma_true),'Copy coefficients changed.');
    assert(norm(t.A_supra_true-t.A_supra_true','fro')<1e-10);
    assert(all(eig(.5*(t.Sigma+t.Sigma'))>-1e-9),'Non-PSD signal covariance.');
end
for split=["validation","test"]
    for trial=1:(cfg.validation_count*(split=="validation")+cfg.test_count*(split=="test"))
        Q=index(lower(index.split)==split&index.trial==trial,:);
        assert(isscalar(unique(Q.innovation_seed)), ...
            'Innovations are not paired across rho.');
    end
end
report=struct('valid',true,'anchor_p',cfg.anchor_p,'rho_count',numel(truths), ...
    'record_count',height(index),'paired_innovations',true);
end
