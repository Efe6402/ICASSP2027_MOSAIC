function [l, err] = PGL_solver(P, q, C, d, tol, rho)
%PGL_SOLVER Output-equivalent callable extraction of PGL_graphs.m (2021).
% The released loop freezes l once the residual is <=tol but continues
% recomputing that unchanged residual until iteration 20,000.  We terminate
% at that point: the returned optimizer is exactly the same, while a large
% hyperparameter study avoids millions of mathematically null iterations.

p = diag(P);
mu = zeros(size(C,1),1);
l_aux = p.^(-1).*(C'*mu - q);
l = max(0, l_aux);

k = 1;
res = norm(C*l - d);
err = nan(20000,1);
err(1) = res;
while k < 20000 && res > tol
    mu = mu - rho*(C*l - d);
    l_aux = p.^(-1).*(C'*mu - q);
    l = max(0, l_aux);
    k = k+1;
    res = norm(C*l - d);
    err(k) = res;
end
err=err(1:k);
end
