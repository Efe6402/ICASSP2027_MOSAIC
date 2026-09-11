%% End-to-end protocol check
clearvars; close all; clc;
root=fileparts(mfilename('fullpath')); addpath(root);
cfg=sample_size_config('quick');
summary=run_model_selection(cfg);
assignin('base','quick_test_summary',summary);
fprintf('\nQuick-test results:\n%s\n',summary.output_dir);
