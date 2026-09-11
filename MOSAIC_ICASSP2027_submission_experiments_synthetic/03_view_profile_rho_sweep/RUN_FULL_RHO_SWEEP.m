%% View-profile mixing sweep
clearvars; close all; clc;
root=fileparts(mfilename('fullpath')); addpath(root);
cfg=rho_sweep_config('full');
summary=run_rho_sweep(cfg);
assignin('base','rho_sweep_summary',summary);
fprintf('\nResults:\n%s\n',summary.output_dir);
