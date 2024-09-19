function [T,P] = right_compression_matrix(n,k)

%RIGHT_COMPRESSION_MATRIX  Right compression matrix for rectangular MEP
%
% [T,P] = RIGHT_COMPRESSION_MATRIX(n,k) returns right compression matrix for the 
% rectangular multiparameter eigenvalue problem and operator determinant method
%
% A{1} x + lambda(1) A{2} x + ... + lambda(k) A{k+1} x = 0 
%
% where matrices are of size m x n and m = n+k-1

% Bor Plestenjak, 2022

% See: M.E.Hochstenbach, T.Kosir, B.Plestenjak: Numerical methods for rectangular 
% multiparameter eigenvalue problems, with applications to finding optimal 
% ARMA and LTI models. Numer Linear Algebra Appl. 2023; e2540

% Based on: F. F. Alsubaie: H2 Optimal Model Reduction for Linear Dynamic 
% Systems and the Solution of Multiparameter Matrix Pencil Problems, PhD,
% Imperial College London, 2019.

P = zeros(n^k,k);
q = (0:n-1)';
z = ones(n,1);

for j = 1:k
    if j==1
        vec = q;
    else
        vec = z;
    end
    for s=2:k
        if j==s
            vec = kron(vec,q);
        else
            vec = kron(vec,z);
        end
    end
    P(:,j) = vec;
end

fac = n.^((k-1:-1:0)');
Q = zeros(n^k,1);
for j = 1:n^k
    Q(j) = sort(P(j,:))*fac;
end

x = unique(Q);
T = spalloc(n^k,length(x),n^k);
for j = 1:length(x)
    vrs = find(Q==x(j));
    T(vrs,j) = 1;
end