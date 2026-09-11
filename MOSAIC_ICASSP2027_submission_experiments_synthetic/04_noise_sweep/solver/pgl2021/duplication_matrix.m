function Dm = duplication_matrix(N)
% function M = duplication_matrix(N)
% M = DuplicationM(n)
% Return duplication matrix order N, with -1 along the off-diagonals of
% the symmetric matrix A
%
%
% Ouput are sparse
%
% duplication_matrix(size(A))*abs(A(itril(size(A)))) == vec(A) 
%
%
% Based on the code from 
% Ref : Magnus, Jan R.; Neudecker, Heinz (1980), "The elimination matrix:
% some lemmas and applications", Society for Industrial and Applied Mathematics.
% Journal on Algebraic and Discrete Methods 1 (4): 422�449,  
% doi:10.1137/0601049, ISSN 0196-5212.
% Author: Sundeep Prabhakar Chepuri <spchepuri@iisc.ac.in>
% Date: 17/Aug/2020

if isscalar(N)
    N = [N N];
end

[I, J] = itril(N);

% Find the sub/sup diagonal part that can flip to other side
loctri = find(I~=J & J<=N(1) & I<=N(2));
% Find the main diagonal part 
locdiag = find(I==J);

% Convert to linear indices
Idiag = sub2ind(N,I(locdiag), J(locdiag));
I1 =  sub2ind(N, I(loctri), J(loctri));
% Indices of the flipped part
Itransposed = sub2ind(N, J(loctri), I(loctri));
       
Dm = sparse([Idiag; I1; Itransposed], ...
           [locdiag;loctri;loctri], ...
           [ones(length(Idiag),1);-ones(length(I1)+length(Itransposed),1)]...
           , prod(N), N(1)*(N(1)+1)/2);
end
