%DEMO_RECT_MULTIPAREIG   demo rectangular multiparameter eigenvalue problem
%
% We solve a three-parameter rectangular eigenvalue problem
%
% (A1 + lambda_1 A2 + lambda_2 A3 + lambda_3 A4) x = 0
%
% where A1, A2, A3, A4 are matrices of size (n+2) x n
%
% Such problem has n*(n+1)*(n+2)/6 eigenvalues
% 
% This is demo for Algorithms 1 and 2 for rectangular multiparameter 
% eigenvalue problems in Section 3 of M.E.Hochstenbach, T.Kosir, 
% B.Plestenjak: On the solution of rectangular multiparameter eigenvalue 
% problems, arXiv 2212.01867

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 6.11.2022

n = 10;
A = cell(1,4);

for k = 1:4
    A{k} = randn(n+2,n);
end

% first approach - using compressed Delta matrices 
tic; [lambda1,X1] = rect_multipareig(A); t1 = toc;

n1 = size(lambda1,1);
err1 = [];
for j = 1:n1
    M = A{1};
    for k = 1:3
        M = M + lambda1(j,k)*A{k+1};
    end
    err1(j,1) = norm(M*X1(:,j));
end

fprintf('First approach, found %d eigenvalues, time: %7.1e, max residual: %7.1e\n', n1,t1,max(err1));

% second approach - using random matrices and conversion to square MEP
opts = [];
opts.method = 'mep';
tic; [lambda2,X2] = rect_multipareig(A,opts); t2 = toc;

n2 = size(lambda2,1);
err2 = [];
for j = 1:n2
    M = A{1};
    for k = 1:3
        M = M + lambda2(j,k)*A{k+1};
    end
    err2(j,1) = norm(M*X2(:,j));
end

fprintf('Second approach, found %d eigenvalues, time: %7.1e, max residual: %7.1e\n', n2,t2,max(err2));
        