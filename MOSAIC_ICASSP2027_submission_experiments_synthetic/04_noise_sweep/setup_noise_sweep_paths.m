function setup_noise_sweep_paths(root)
%SETUP_NOISE_SWEEP_PATHS Add only package-local solver paths.
if nargin<1||isempty(root), root=fileparts(mfilename('fullpath')); end
addpath(root);
addpath(fullfile(root,'solver','mosaic'));
addpath(fullfile(root,'solver','pgl2021'));
addpath(fullfile(root,'solver','zw2024'));
end
