function [truth,index,levels]=prepare_noise_signal_bank(cfg)
%PREPARE_NOISE_SIGNAL_BANK Build paired view-wise kappa x SNR observations.
%
% A single clean realization Yclean and a single standard-Gaussian innovation
% E0 are drawn per trial.  The SAME E0 is reused for every kappa and SNR.
% For view k, the prescribed variance multiplier is
%
%   v_k(kappa) = 1 + kappa * (v_k^* - 1),
%
% where v^* is cfg.view_variance_endpoint.  Thus kappa=0 is the equal-
% variance IID/AWGN endpoint and kappa=1 is the original heterogeneous-view
% endpoint.  No node-dependent gains are used.
if ~isfolder(cfg.data_dir), mkdir(cfg.data_dir); end
indexfile=fullfile(cfg.data_dir,'BANK_INDEX.csv');
levelsfile=fullfile(cfg.data_dir,'NOISE_CONDITIONS.csv');
expected=numel(cfg.snr_db)*numel(cfg.kappa_grid)*cfg.test_count;

if ~cfg.overwrite_data&&isfile(indexfile)&&isfile(levelsfile)
    index=readtable(indexfile,'TextType','string');
    levels=readtable(levelsfile,'TextType','string');
    good=height(index)==expected&& ...
        height(levels)==numel(cfg.snr_db)*numel(cfg.kappa_grid)&& ...
        all(index.sample_count==cfg.anchor_p)&&all(isfile(index.file));
    if good
        expected_kappa=repelem(cfg.kappa_grid(:),numel(cfg.snr_db),1);
        expected_snr=repmat(cfg.snr_db(:),numel(cfg.kappa_grid),1);
        good=good&&max(abs(double(levels.kappa)-expected_kappa))<1e-12&& ...
            all((isinf(double(levels.target_snr_db))&isinf(expected_snr)) | ...
            abs(double(levels.target_snr_db)-expected_snr)<1e-12);
    end
    if good
        endpoint=[levels.endpoint_view_variance_1(1), ...
            levels.endpoint_view_variance_2(1), ...
            levels.endpoint_view_variance_3(1), ...
            levels.endpoint_view_variance_4(1)];
        good=max(abs(endpoint-cfg.view_variance_endpoint))<1e-12;
    end
    if good, truth=load_default_truth(cfg); return; end
end

truth=load_default_truth(cfg);
[V,D]=eig(.5*(truth.Sigma+truth.Sigma'));
F=real(V*diag(sqrt(max(real(diag(D)),0))));

% Cartesian condition table: each kappa is evaluated at every SNR.
nKappa=numel(cfg.kappa_grid); nSnr=numel(cfg.snr_db);
kappa_index=repelem((1:nKappa)',nSnr,1);
snr_index=repmat((1:nSnr)',nKappa,1);
kappa=repelem(cfg.kappa_grid(:),nSnr,1);
target_snr_db=repmat(cfg.snr_db(:),nKappa,1);
condition=(1:numel(kappa))';
noise_amplitude_ratio=zeros(size(target_snr_db));
finite_level=isfinite(target_snr_db);
noise_amplitude_ratio(finite_level)=10.^(-target_snr_db(finite_level)/20);
noise_to_signal_power_ratio=noise_amplitude_ratio.^2;
noise_energy_fraction=noise_to_signal_power_ratio./(1+noise_to_signal_power_ratio);

P=zeros(numel(kappa),truth.K);
for c=1:numel(kappa)
    P(c,:)=view_variance_profile_from_kappa(kappa(c),cfg.view_variance_endpoint);
end
levels=table(condition,kappa_index,kappa,snr_index,target_snr_db, ...
    noise_amplitude_ratio,noise_to_signal_power_ratio,noise_energy_fraction, ...
    repmat("viewwise_kappa_heteroscedastic_gaussian",numel(kappa),1), ...
    P(:,1),P(:,2),P(:,3),P(:,4), ...
    repmat(cfg.view_variance_endpoint(1),numel(kappa),1), ...
    repmat(cfg.view_variance_endpoint(2),numel(kappa),1), ...
    repmat(cfg.view_variance_endpoint(3),numel(kappa),1), ...
    repmat(cfg.view_variance_endpoint(4),numel(kappa),1), ...
    abs(kappa)<1e-12,abs(kappa-1)<1e-12,isinf(target_snr_db),target_snr_db==0, ...
    'VariableNames',{'level','kappa_index','kappa','snr_index','target_snr_db', ...
    'noise_amplitude_ratio','noise_to_signal_power_ratio','noise_energy_fraction', ...
    'noise_model','unrotated_view_variance_1','unrotated_view_variance_2', ...
    'unrotated_view_variance_3','unrotated_view_variance_4', ...
    'endpoint_view_variance_1','endpoint_view_variance_2', ...
    'endpoint_view_variance_3','endpoint_view_variance_4', ...
    'is_iid_kappa_endpoint','is_full_heteroscedastic_endpoint', ...
    'is_clean_endpoint','is_equal_power_endpoint'});

rows=table();
for trial=1:cfg.test_count
    clean_seed=cfg.seed_base+200000+trial;
    noise_seed=cfg.seed_base+50000000+trial;
    old=rng; cleaner=onCleanup(@()rng(old));
    rng(clean_seed,'twister'); Yclean=F*randn(truth.N,cfg.anchor_p);
    rng(noise_seed,'twister'); E0=randn(size(Yclean));
    clean_power=mean(Yclean.^2,'all');
    assert(clean_power>0);

    if cfg.rotate_variance_profile
        profile_rotation=mod(trial-1,truth.K);
    else
        profile_rotation=0;
    end

    for ki=1:nKappa
        kap=cfg.kappa_grid(ki);
        unrotated=view_variance_profile_from_kappa(kap,cfg.view_variance_endpoint);
        applied=circshift(unrotated,[0 profile_rotation]);
        row_variances=kron(applied(:),ones(truth.n,1));

        % View-wise only: every physical node inside a given view receives
        % the same prescribed variance multiplier.  There is no node-wise
        % low/high-noise partition anywhere in this package.
        E=bsxfun(@times,E0,sqrt(row_variances));
        raw_noise_power=mean(E.^2,'all');
        assert(raw_noise_power>0);

        % Realization-wise calibration: preserve the relative view profile
        % while setting the overall pre-SNR noise power equal to Yclean.
        E=E*sqrt(clean_power/raw_noise_power);
        base_view_noise_power=zeros(1,truth.K);
        for vk=1:truth.K
            Ik=(vk-1)*truth.n+(1:truth.n);
            base_view_noise_power(vk)=mean(E(Ik,:).^2,'all');
        end

        for si=1:nSnr
            lev=(ki-1)*nSnr+si;
            alpha=levels.noise_amplitude_ratio(lev);
            Ys=Yclean; Yn=alpha*E; Y=Ys+Yn;
            signal_power=mean(Ys.^2,'all');
            noise_power=mean(Yn.^2,'all');
            observed_power=mean(Y.^2,'all');
            if noise_power==0
                realized_snr=Inf;
            else
                realized_snr=10*log10(signal_power/noise_power);
            end
            realized_power_ratio=noise_power/signal_power;
            realized_fraction=noise_power/(signal_power+noise_power);

            meta=struct('split','test','level',lev,'trial',trial, ...
                'seed',clean_seed,'clean_seed',clean_seed,'noise_seed',noise_seed, ...
                'noise_model','viewwise_kappa_heteroscedastic_gaussian', ...
                'kappa_index',ki,'kappa',kap,'snr_index',si, ...
                'target_snr_db',levels.target_snr_db(lev), ...
                'noise_amplitude_ratio',alpha, ...
                'noise_to_signal_power_ratio',levels.noise_to_signal_power_ratio(lev), ...
                'unrotated_view_variances',unrotated, ...
                'applied_view_variances',applied, ...
                'variance_profile_rotation',profile_rotation, ...
                'base_view_noise_power',base_view_noise_power, ...
                'realized_snr_db',realized_snr, ...
                'realized_noise_to_signal_power_ratio',realized_power_ratio, ...
                'realized_noise_energy_fraction',realized_fraction, ...
                'clean_power',clean_power,'signal_component_power',signal_power, ...
                'noise_component_power',noise_power,'observed_power',observed_power);
            record=make_signal_record(Y,truth,meta);

            folder=fullfile(cfg.data_dir,kappa_folder(kap),snr_folder(levels.target_snr_db(lev)));
            if ~isfolder(folder), mkdir(folder); end
            file=fullfile(folder,sprintf('test_trial%02d.mat',trial));
            save(file,'record','-v7.3');

            row=table(lev,ki,kap,si,levels.target_snr_db(lev),alpha, ...
                levels.noise_to_signal_power_ratio(lev),levels.noise_energy_fraction(lev), ...
                cfg.anchor_p,trial,clean_seed,noise_seed,profile_rotation, ...
                unrotated(1),unrotated(2),unrotated(3),unrotated(4), ...
                applied(1),applied(2),applied(3),applied(4), ...
                alpha^2*base_view_noise_power(1),alpha^2*base_view_noise_power(2), ...
                alpha^2*base_view_noise_power(3),alpha^2*base_view_noise_power(4), ...
                clean_power,signal_power,noise_power,observed_power, ...
                realized_power_ratio,realized_fraction,realized_snr,string(file), ...
                'VariableNames',{'level','kappa_index','kappa','snr_index', ...
                'target_snr_db','noise_amplitude_ratio','noise_to_signal_power_ratio', ...
                'noise_energy_fraction','sample_count','trial','clean_seed','noise_seed', ...
                'variance_profile_rotation','unrotated_view_variance_1', ...
                'unrotated_view_variance_2','unrotated_view_variance_3', ...
                'unrotated_view_variance_4','applied_view_variance_1', ...
                'applied_view_variance_2','applied_view_variance_3', ...
                'applied_view_variance_4','realized_view_noise_power_1', ...
                'realized_view_noise_power_2','realized_view_noise_power_3', ...
                'realized_view_noise_power_4','clean_power','signal_component_power', ...
                'noise_component_power','observed_power', ...
                'realized_noise_to_signal_power_ratio','realized_noise_energy_fraction', ...
                'realized_snr_db','file'});
            rows=append_rows(rows,row);
        end
    end
    clear cleaner Yclean E0 E
end
index=sortrows(rows,{'level','trial'});
writetable(index,indexfile);
writetable(levels,levelsfile);
end

function s=kappa_folder(kappa)
s=sprintf('kappa_%03d',round(100*kappa));
end

function s=snr_folder(snr_db)
if isinf(snr_db), s='snr_clean'; else, s=sprintf('snr_%02d',round(snr_db)); end
end

function T=append_rows(T,S)
if isempty(T), T=S; else, T=[T;S]; end
end
