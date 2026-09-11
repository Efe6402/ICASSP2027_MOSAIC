function [C] = KronSum(A,B)
% C = kron(A,eye(size(B))) + kron(eye(size(A)),B);
C = kron(eye(size(B)) , A) + kron(B, eye(size(A)));
end