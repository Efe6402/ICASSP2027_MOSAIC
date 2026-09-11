%% p sweep with p=5000-selected configurations and validation thresholds
clearvars; close all; clc;
root=fileparts(mfilename('fullpath')); addpath(root);
cfg=sample_size_config('full');
summary=run_p5000_fixed_transfer(cfg);
assignin('base','p5000_transfer_summary',summary);
fprintf('\nResults:\n%s\n',summary.output_dir);
