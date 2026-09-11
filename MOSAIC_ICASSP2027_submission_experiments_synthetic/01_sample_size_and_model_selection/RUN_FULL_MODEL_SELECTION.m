%% Validation model selection at every sample count
clearvars; close all; clc;
root=fileparts(mfilename('fullpath')); addpath(root);
cfg=sample_size_config('full');
summary=run_model_selection(cfg);
assignin('base','model_selection_summary',summary);
fprintf('\nResults:\n%s\n',summary.output_dir);
