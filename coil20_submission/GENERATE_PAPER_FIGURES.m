%% Regenerate the four paper figures from the packaged selected model
clearvars; clc;

package_root = fileparts(mfilename('fullpath'));
fit_file = fullfile(package_root,'selected_top1', ...
    'top01_r3_profile2508_G08_beta_c_010.mat');
data_file = fullfile(package_root,'data','global_6view_phases', ...
    'coil20_global6_phase_00.mat');
output_dir = fullfile(package_root,'paper_figures');

generate_coil20_paper_figures(fit_file,data_file,output_dir);
fprintf('\nPaper figures saved to: %s\n',output_dir);
