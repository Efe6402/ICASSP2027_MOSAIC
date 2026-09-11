function report=validate_noise_bank(cfg,truth,index,levels)
%VALIDATE_NOISE_BANK Verify view-wise kappa interpolation, pairing, and SNR.
nKappa=numel(cfg.kappa_grid); nSnr=numel(cfg.snr_db);
expected=nKappa*nSnr*cfg.test_count;
assert(height(index)==expected&&height(levels)==nKappa*nSnr);
assert(all(index.sample_count==cfg.anchor_p)&&all(isfile(index.file)));

expected_kappa=repelem(cfg.kappa_grid(:),nSnr,1);
expected_snr=repmat(cfg.snr_db(:),nKappa,1);
assert(max(abs(double(levels.kappa)-expected_kappa))<1e-12);
assert(all((isinf(double(levels.target_snr_db))&isinf(expected_snr)) | ...
    abs(double(levels.target_snr_db)-expected_snr)<1e-12));
endpoint=[levels.endpoint_view_variance_1(1),levels.endpoint_view_variance_2(1), ...
    levels.endpoint_view_variance_3(1),levels.endpoint_view_variance_4(1)];
assert(max(abs(endpoint-cfg.view_variance_endpoint))<1e-12);
assert(abs(mean(endpoint)-1)<1e-12&&std(endpoint)>0);

% The prescribed unrotated variance profile must exactly follow the kappa law.
for ki=1:nKappa
    kap=cfg.kappa_grid(ki);
    Q=levels(levels.kappa_index==ki,:);
    expected_profile=view_variance_profile_from_kappa(kap,cfg.view_variance_endpoint);
    got=[Q.unrotated_view_variance_1(1),Q.unrotated_view_variance_2(1), ...
        Q.unrotated_view_variance_3(1),Q.unrotated_view_variance_4(1)];
    assert(max(abs(got-expected_profile))<1e-12);
    assert(all(Q.unrotated_view_variance_1==got(1))&& ...
        all(Q.unrotated_view_variance_2==got(2))&& ...
        all(Q.unrotated_view_variance_3==got(3))&& ...
        all(Q.unrotated_view_variance_4==got(4)));
end
assert(max(abs(view_variance_profile_from_kappa(0,cfg.view_variance_endpoint)-ones(1,4)))<1e-12);
assert(max(abs(view_variance_profile_from_kappa(1,cfg.view_variance_endpoint)-cfg.view_variance_endpoint))<1e-12);

assert(truth.N==truth.n*truth.K&&isequal(size(truth.Sigma),[truth.N truth.N]));
assert(norm(truth.A_supra_true-truth.A_supra_true','fro')<1e-10);
assert(all(eig(.5*(truth.Sigma+truth.Sigma'))>-1e-9));

for trial=1:cfg.test_count
    Qt=index(index.trial==trial,:);
    assert(height(Qt)==nKappa*nSnr);
    assert(isscalar(unique(Qt.clean_seed))&&isscalar(unique(Qt.noise_seed)), ...
        'Each trial must use one clean realization and one base Gaussian innovation across the full grid.');
    assert(isscalar(unique(Qt.variance_profile_rotation)));

    for ki=1:nKappa
        kap=cfg.kappa_grid(ki);
        Q=Qt(Qt.kappa_index==ki,:);
        assert(height(Q)==nSnr);
        unrotated=view_variance_profile_from_kappa(kap,cfg.view_variance_endpoint);
        if cfg.rotate_variance_profile
            expected_profile=circshift(unrotated,[0 mod(trial-1,truth.K)]);
        else
            expected_profile=unrotated;
        end
        applied=[Q.applied_view_variance_1(1),Q.applied_view_variance_2(1), ...
            Q.applied_view_variance_3(1),Q.applied_view_variance_4(1)];
        assert(max(abs(applied-expected_profile))<1e-12);
        assert(all(Q.applied_view_variance_1==applied(1))&& ...
            all(Q.applied_view_variance_2==applied(2))&& ...
            all(Q.applied_view_variance_3==applied(3))&& ...
            all(Q.applied_view_variance_4==applied(4)), ...
            'Applied VIEW variance profile must stay fixed across SNR for a given trial and kappa.');

        scale=max(1,max(abs(Q.clean_power)));
        assert(max(abs(Q.signal_component_power-Q.clean_power))<1e-12*scale, ...
            'Direct additive noise must not attenuate the clean signal.');
        assert(max(abs(Q.realized_noise_to_signal_power_ratio- ...
            Q.noise_to_signal_power_ratio))<1e-12);
        assert(max(abs(Q.realized_noise_energy_fraction-Q.noise_energy_fraction))<1e-12);
        finite=isfinite(Q.target_snr_db);
        assert(max(abs(Q.realized_snr_db(finite)-Q.target_snr_db(finite)))<1e-10);
    end
end

assert(all(index.noise_component_power(isinf(index.target_snr_db))==0));
assert(all(abs(index.noise_to_signal_power_ratio(index.target_snr_db==0)-1)<1e-12));

% Explicit endpoint checks: kappa=0 has equal prescribed view variances;
% kappa=1 has the heterogeneous endpoint (up to trial rotation).
Q0=index(abs(index.kappa)<1e-12,:);
applied0=double(Q0{:,{'applied_view_variance_1','applied_view_variance_2', ...
    'applied_view_variance_3','applied_view_variance_4'}});
assert(max(abs(applied0(:)-1))<1e-12,'kappa=0 must be the equal-variance IID endpoint.');
Q1=index(abs(index.kappa-1)<1e-12&index.target_snr_db==0,:);
view_powers=double(Q1{:,{'realized_view_noise_power_1','realized_view_noise_power_2', ...
    'realized_view_noise_power_3','realized_view_noise_power_4'}});
assert(all(max(view_powers,[],2)>min(view_powers,[],2)), ...
    'kappa=1 must produce unequal view-wise noise powers.');

z=load(char(index.file(1)),'record');
assert(z.record.n==truth.n&&z.record.K==truth.K&&z.record.sample_count==cfg.anchor_p);
report=struct('valid',true,'sample_count',cfg.anchor_p, ...
    'kappa_count',nKappa,'snr_count',nSnr,'condition_count',height(levels), ...
    'test_record_count',height(index),'paired_clean_signal',true, ...
    'paired_base_gaussian_innovation_across_kappa_and_snr',true, ...
    'direct_additive_noise',true,'view_dependent_only',true, ...
    'node_dependent_noise',false,'kappa_zero_equal_variance',true, ...
    'kappa_one_endpoint_profile',true, ...
    'endpoint_view_noise_variance_1',cfg.view_variance_endpoint(1), ...
    'endpoint_view_noise_variance_2',cfg.view_variance_endpoint(2), ...
    'endpoint_view_noise_variance_3',cfg.view_variance_endpoint(3), ...
    'endpoint_view_noise_variance_4',cfg.view_variance_endpoint(4), ...
    'rotated_profile_across_trials',cfg.rotate_variance_profile, ...
    'minimum_finite_snr_db',min(cfg.snr_db(isfinite(cfg.snr_db))));
end
