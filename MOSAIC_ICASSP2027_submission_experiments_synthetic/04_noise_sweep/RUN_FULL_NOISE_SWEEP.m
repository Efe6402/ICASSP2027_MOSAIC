%% Full fixed-configuration VIEW-wise kappa x SNR Gaussian-noise experiment
clearvars; clc;
root=fileparts(mfilename('fullpath')); addpath(root);
cfg=noise_sweep_config('full');
summary=run_noise_sweep(cfg);
assignin('base','noise_sweep_summary',summary);
fprintf('\nKappa x SNR noise-sweep results:\n%s\n',summary.output_dir);
