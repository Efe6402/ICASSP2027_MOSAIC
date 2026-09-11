%% Model selection followed by p=5000 fixed-configuration transfer
clearvars; close all; clc;
root=fileparts(mfilename('fullpath')); addpath(root);
cfg=sample_size_config('full');
selection_summary=run_model_selection(cfg);
selected_dir=fullfile(selection_summary.output_dir,'p_05000');
transfer_cfg=sample_size_config('full',struct('selected_config_dir',selected_dir));
transfer_summary=run_p5000_fixed_transfer(transfer_cfg);
assignin('base','model_selection_summary',selection_summary);
assignin('base','p5000_transfer_summary',transfer_summary);
fprintf('\nModel selection:\n%s\n',selection_summary.output_dir);
fprintf('p=5000 transfer:\n%s\n',transfer_summary.output_dir);
