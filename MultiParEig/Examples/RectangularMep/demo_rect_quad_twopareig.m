%DEMO_RECT_QUAD_TWOPAREIG   demo rectangular quadratic twoparameter eigenvalue problem
%
% We solve a rectangular quadratic two-parameter eigenvalue problem
%
% (A + lambda*B + mu*C + lambda^2*D+ lambda*mu*E + mu^2*F)x = 0
% 
% where A, B, C, D, E, F are (n+1) x n matrices
%
% Such problem has 2*n*(n+1) eigenvalues
%
% This is demo for algorithms for rectangular quadratic two-parameter 
% eigenvalue problem in Section 5 of M.E.Hochstenbach, T.Kosir, 
% B.Plestenjak: On the solution of rectangular multiparameter eigenvalue 
% problems, arXiv 2212.01867

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 6.11.2022

n = 10;

A = randn(n+1,n);
B = randn(n+1,n);
C = randn(n+1,n);
D = randn(n+1,n);
E = randn(n+1,n);
F = randn(n+1,n);

% first approach - using compressed Delta matrices 
opts = [];
opts.method = 'mep';
opts.showrank  = 1;
tic; [lambda1,mu1,X1] = rect_quad_twopareig(A,B,C,D,E,F,opts); t1 = toc;

n1 = size(lambda1,1);
err1 = [];
for j = 1:n1
    M = A + lambda1(j)*B +mu1(j)*C + lambda1(j)^2*D + lambda1(j)*mu1(j)*E + mu1(j)^2*F;
    err1(j,1) = norm(M*X1(:,j));
end

fprintf('First approach, found %d eigenvalues, time: %7.1e, max residual: %7.1e\n', n1,t1,max(err1));

% second approach - linearize to linear RMEP
opts = [];
opts.method = 'linearize';
opts.showrank  = 1;
tic; [lambda2,mu2,X2] = rect_quad_twopareig(A,B,C,D,E,F,opts); t2 = toc;

n2 = size(lambda2,1);
err2 = [];
for j = 1:n2
    M = A + lambda2(j)*B +mu2(j)*C + lambda2(j)^2*D + lambda2(j)*mu2(j)*E + mu2(j)^2*F;
    err2(j,1) = norm(M*X2(:,j));
end

fprintf('Second approach, found %d eigenvalues, time: %7.1e, max residual: %7.1e\n', n2,t2,max(err2));

% third approach - linearize to linear RMEP and Vandermonde compression
opts = [];
opts.method = 'compress';
opts.showrank  = 1;
tic; [lambda3,mu3,X3] = rect_quad_twopareig(A,B,C,D,E,F,opts); t3 = toc;

n3 = size(lambda3,1);
err3 = [];
for j = 1:n3
    M = A + lambda3(j)*B +mu3(j)*C + lambda3(j)^2*D + lambda3(j)*mu3(j)*E + mu3(j)^2*F;
    err3(j,1) = norm(M*X3(:,j));
end

fprintf('Third approach, found %d eigenvalues, time: %7.1e, max residual: %7.1e\n', n3,t3,max(err3));

        