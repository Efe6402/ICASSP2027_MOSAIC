%% MOSAIC S&P 500 experiment
clearvars; close all; clc;

root = fileparts(mfilename('fullpath'));
addpath(root);

stamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
workRoot = fullfile(root,'work',['run_' stamp]);
resultRoot = fullfile(root,'results',['run_' stamp]);

cfg = struct();
cfg.output_root = workRoot;
cfg.stocks_dir = fullfile(root,'data','stocks');
cfg.prepared_target_mat = fullfile(root,'data','prepared_offset0', ...
    'mosaic_finance_multihorizon_dataset.mat');
cfg.resume = false;

run_sp500_submission_core(cfg);
curate_sp500_submission_outputs(workRoot,resultRoot);

if isfolder(workRoot)
    rmdir(workRoot,'s');
end

fprintf('\nSubmission outputs:\n%s\n',resultRoot);
