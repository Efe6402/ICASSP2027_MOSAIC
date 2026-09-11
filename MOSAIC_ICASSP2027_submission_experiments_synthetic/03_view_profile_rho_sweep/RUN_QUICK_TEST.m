%% End-to-end rho protocol check
clearvars; close all; clc;
root=fileparts(mfilename('fullpath')); addpath(root);
cfg=rho_sweep_config('quick');
summary=run_rho_sweep(cfg);
assignin('base','rho_quick_summary',summary);
fprintf('\nQuick-test results:\n%s\n',summary.output_dir);
