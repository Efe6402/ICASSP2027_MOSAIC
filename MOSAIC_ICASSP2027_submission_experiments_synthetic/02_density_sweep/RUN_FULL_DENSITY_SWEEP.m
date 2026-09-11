%% Exact-count modal-density sweep
clearvars; close all; clc;
root=fileparts(mfilename('fullpath')); addpath(root);
cfg=density_sweep_config('full');
summary=run_density_sweep(cfg);
assignin('base','density_sweep_summary',summary);
fprintf('\nResults:\n%s\n',summary.output_dir);
