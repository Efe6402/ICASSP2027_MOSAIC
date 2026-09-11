function v=view_variance_profile_from_kappa(kappa,endpoint)
%VIEW_VARIANCE_PROFILE_FROM_KAPPA Interpolate IID -> heterogeneous views.
%
%   v(kappa) = 1 + kappa * (endpoint - 1)
%
% The endpoint is required to have unit mean, so every interpolated profile
% also has unit mean.  These are VARIANCE multipliers; samples are scaled by
% sqrt(v) when the Gaussian innovation is constructed.
assert(isscalar(kappa)&&isfinite(kappa)&&kappa>=0&&kappa<=1);
endpoint=double(endpoint(:)');
assert(numel(endpoint)==4&&all(endpoint>0)&&abs(mean(endpoint)-1)<1e-12);
v=ones(size(endpoint))+double(kappa)*(endpoint-1);
assert(all(v>0)&&abs(mean(v)-1)<1e-12);
end
