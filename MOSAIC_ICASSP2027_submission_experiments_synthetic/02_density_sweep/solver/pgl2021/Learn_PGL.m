function [Lp_i, Lq_i] = Learn_PGL(Sp, Sq, param)
%LEARN_PGL Callable extraction of the sparse-factor solver in PGL_graphs.m.
% Algebra, constants, initialization, and 1e-4 output threshold are kept
% identical to the public 2021 MATLAB implementation.

p = size(Sp,1);
q = size(Sq,1);

Dp = duplication_matrix(p);
Dq = duplication_matrix(q);

P0 = blkdiag(2*param.b1*Dp'*Dp, 2*param.b2*Dq'*Dq);
q1p = [vec(Sp)'*Dp]';
q1q = [vec(Sq)'*Dq]';
q0 = [q1p; q1q];

Cp = [vec(eye(p))'*Dp; kron(ones(p,1)',eye(p))*Dp];
dp = [p; zeros(p,1)];
Cq = [vec(eye(q))'*Dq; kron(ones(q,1)',eye(q))*Dq];
dq = [q; zeros(q,1)];
C = blkdiag(Cp, Cq);
d = [dp; dq];

[l_wf, ~] = PGL_solver(P0, q0, C, d, 1e-6, 0.0051);
v1 = size(Dp,2);
Lp_i = full(reshape(Dp*l_wf(1:v1),p,p));
Lp_i(abs(Lp_i)<1e-4) = 0;
Lq_i = full(reshape(Dq*l_wf(v1+1:end),q,q));
Lq_i(abs(Lq_i)<1e-4) = 0;
end
