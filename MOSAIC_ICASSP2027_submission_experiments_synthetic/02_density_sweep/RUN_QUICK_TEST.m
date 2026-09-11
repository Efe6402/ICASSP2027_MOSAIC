%% End-to-end density protocol check
clearvars; close all; clc;
root=fileparts(mfilename('fullpath')); addpath(root);
cfg=density_sweep_config('quick');
summary=run_density_sweep(cfg);
assignin('base','density_quick_summary',summary);
fprintf('\nQuick-test results:\n%s\n',summary.output_dir);
