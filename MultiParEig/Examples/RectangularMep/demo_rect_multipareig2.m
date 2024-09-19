%DEMO_RECT_MULTIPAREIG2   demo rectangular multiparameter eigenvalue problem
%
% We solve a multiparameter rectangular eigenvalue problem
%
% (A1 + lambda_1 A2 + ... + lambda_k A_{p+1}) x = 0
%
% where A1, A2, A3, A_{p+1} are matrices of size (n+p-1) x n
%
% Such problem has nchoosek(n+p-1,p) eigenvalues
% 
% This is demo for rectangular multiparameter eigenvalue problems in
% M.E.Hochstenbach, T.Kosir, B.Plestenjak: Numerical methods for rectangular 
% multiparameter eigenvalue problems, with applications to finding optimal 
% ARMA and LTI models. Numer Linear Algebra Appl. 2023; e2540

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 28.08.2024

p = 6;  % number of parameters
n = 5;  % size of matrices is (n+p-1) x n
A = cell(1,p+1);

for k = 1:p+1
    A{k} = randn(n+p-1,n);
end

tic; [lambda1,X1] = rect_multipareig(A); t1 = toc;

n1 = size(lambda1,1);
err1 = [];
for j = 1:n1
    M = A{1};
    for k = 1:p
        M = M + lambda1(j,k)*A{k+1};
    end
    err1(j,1) = norm(M*X1(:,j));
end

fprintf('Found %d eigenvalues, time: %7.1e, max residual: %7.1e\n', n1,t1,max(err1));
