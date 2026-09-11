function S = compute_zhang_wai_2024_distance_matrix(Y)
%COMPUTE_ZHANG_WAI_2024_DISTANCE_MATRIX Paper Eq. (10), exactly uncentered.
% S_ij = (1/M) sum_m |y_i^(m)-y_j^(m)|^2.

M=size(Y,2);
if M<1, error('Y must contain at least one graph-signal observation.'); end
q=sum(Y.^2,2)/M;
S=q+q'-2*(Y*Y')/M;
S=0.5*(S+S');
% Only remove roundoff below zero; no rescaling/normalization is performed.
S(S<0 & S>-1e-10)=0;
if any(S(:)<-1e-8), error('Distance matrix has materially negative entries.'); end
S(1:size(S,1)+1:end)=0;
end
