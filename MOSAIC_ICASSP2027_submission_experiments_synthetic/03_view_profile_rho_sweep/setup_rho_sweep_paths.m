function setup_rho_sweep_paths(root)
%SETUP_RHO_SWEEP_PATHS Add this experiment and its three solver folders.
if nargin<1, root=fileparts(mfilename('fullpath')); end
addpath(root,fullfile(root,'solver','mosaic'), ...
    fullfile(root,'solver','pgl2021'),fullfile(root,'solver','zw2024'));
end
