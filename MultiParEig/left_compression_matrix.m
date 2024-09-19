function [T,RWS] = left_compression_matrix(m,k)

%LEFT_COMPRESSION_MATRIX  Left compression matrix for rectangular MEP
%
% [T,RWS] = LEFT_COMPRESSION_MATRIX(m,k) returns left compression matrix T 
% and list of row indices for the rectangular multiparameter eigenvalue 
% problem and operator determinant method
%
% A{1}*x + lambda(1)*A{2}*x + ... + lambda(k)*A{k+1}*x = 0 
%
% where matrices are of size m x n and m = n + k -1
%
% RWS is a list of strictly ordered k-indices with elements from 1:m

% Bor Plestenjak, 2022
% BP 08.11.2023 Faster method

% See: M.E.Hochstenbach, T.Kosir, B.Plestenjak: Numerical methods for rectangular 
% multiparameter eigenvalue problems, with applications to finding optimal 
% ARMA and LTI models. Numer Linear Algebra Appl. 2023; e2540

% Based on: F. F. Alsubaie: H2 Optimal Model Reduction for Linear Dynamic 
% Systems and the Solution of Multiparameter Matrix Pencil Problems, PhD,
% Imperial College London, 2019.

RWS = strictly_ordered(1,m,k);
NN = size(RWS,1);

A = zeros(NN,3);
for row = 1:NN
    col = RWS(row,1)-1;
    for q = 2:k
        col = col*m + RWS(row,q)-1;
    end
    A(row,:) = [row col+1 1];
end

T = sparse(A(:,1),A(:,2),A(:,3),NN,m^k,NN);

function w = strictly_ordered(a,b,k)
if k==1
    w = (a:b)';
elseif b-a < k-1
    w = [];
elseif b-a == k-1
    w = a:b;
else
    w = [];
    for j = a:b-k+1
        subpart = strictly_ordered(j+1,b,k-1);
        w = [w; j*ones(size(subpart,1),1) subpart];
    end
end
