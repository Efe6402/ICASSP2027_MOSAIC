function setup_three_method_paths(root)
%SETUP_THREE_METHOD_PATHS Add this experiment and its three solver folders.
if nargin<1, root=fileparts(mfilename('fullpath')); end
addpath(root,fullfile(root,'solver','mosaic'), ...
    fullfile(root,'solver','pgl2021'),fullfile(root,'solver','zw2024'));
end
