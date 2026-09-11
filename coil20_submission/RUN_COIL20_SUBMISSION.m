%% COIL-20 six-view MOSAIC experiment
clearvars; clc;

package_root = fileparts(mfilename('fullpath'));
addpath(package_root,fullfile(package_root,'solver'));

summary = run_coil20_experiment_core(package_root);
assignin('base','coil20_submission_summary',summary);

fprintf('\nCOIL-20 experiment complete.\nResults: %s\n',summary.output_dir);
