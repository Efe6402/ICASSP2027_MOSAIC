%% Small VIEW-wise kappa x SNR Gaussian-noise end-to-end package check
clearvars; clc;
root=fileparts(mfilename('fullpath')); addpath(root);
cfg=noise_sweep_config('quick',struct('overwrite_data',true));
summary=run_noise_sweep(cfg);
assignin('base','noise_quick_summary',summary);
fprintf('\nQuick-test results:\n%s\n',summary.output_dir);
