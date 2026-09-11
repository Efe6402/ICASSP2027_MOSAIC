function [X1 T S] = DataGen(Lp, Lq, N)

L = KronSum(Lp, Lq);
L = L/trace(L)*size(L,1); 

p = [size(Lp,1) size(Lq,1)];
[V,D] = eig(full(L));
sigma = pinv(D);
mu = zeros(1,size(L,1));
gftcoeff = mvnrnd(mu,sigma, N);
X1 = V*gftcoeff';
T = 0;
S = 0;
for i = 1:N
    x = X1(:,i);
    Xi = reshape(x, p(1), p(2));
    T = T + cov(Xi');
    S = S + cov(Xi);
    
end
end