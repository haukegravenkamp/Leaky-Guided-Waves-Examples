%DEMO_RECT_QUAD_MULTIPAREIG   demo rectangular quadratic multiparameter eigenvalue problem
%
% We solve a rectangular quadratic multiparameter eigenvalue problem for 
% 2, 3, and 4 parameters
% 
% For k=2 (A{1},...,A{6} are (n+1) x n matrices) the form is:
%
% (A{1} + mu(1)*A{2} + mu(2)*A{3} + mu(1)^2*A{4} + mu(1)*mu(2)*A{5} + mu(2)^2*A{6})x = 0
% 
% For k=3 (A{1},...A{10} are (n+2) x n matrices) the form is:
%
% (A{1} + mu(1)*A{2} + mu(2)*A{3} + mu(3)*A{4} + mu(1)^2*A{5} + mu(1)*mu(2)*A{6} 
%       + mu(1)*mu(3)*A{7} + mu(2)^2*A{8} + mu(2)*mu(3)*A{9} + mu(3)^2*A{10})x = 0
% 
% For k=4 (A{1},...A{15} are (n+3) x n matrices) the form is:
%
% (A{1} + mu(1)*A{2} + mu(2)*A{3} + mu(3)*A{4} + mu(4)*A{5} + mu(1)^2*A{6} + mu(1)*mu(2)*A{7} 
%       + mu(1)*mu(3)*A{8} + mu(1)*mu(4)*A{9} + mu(2)^2*A{10} + mu(2)*mu(3)*A{11} 
%       + mu(2)*mu(4)*A{12} + mu(3)^2*A{13} + mu(3)*mu(4)*A{14} + mu(4)^2*A{15})x = 0

% Bor Plestenjak 2024

% See: M.E.Hochstenbach, T.Kosir, B.Plestenjak: Numerical methods for rectangular 
% multiparameter eigenvalue problems, with applications to finding optimal 
% ARMA and LTI models. Numer Linear Algebra Appl. 2023; e2540

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Two parameters, k = 2
n = 10;
m = 6; 
A = cell(m,1);
for j = 1:m
    A{j} = rand(n+1,n);
end

tic; [lambda1,X1] = rect_quad_multipareig(A); t1 = toc;

n1 = size(lambda1,1);
err1 = [];
for j = 1:n1
    mu(1) = lambda1(j,1);
    mu(2) = lambda1(j,2);
    M = A{1} + mu(1)*A{2} + mu(2)*A{3} + mu(1)^2*A{4} + mu(1)*mu(2)*A{5} + mu(2)^2*A{6};
    err1(j,1) = norm(M*X1(:,j));
end

fprintf('k=2, matrices (%d x %d), found %d eigenvalues, time: %7.1e, max residual: %7.1e\n',n+1,n,n1,t1,max(err1));

% Three parameters, k = 3
n = 5;
m = 10; 
A = cell(m,1);
for j = 1:m
    A{j} = rand(n+2,n);
end

tic; [lambda2,X2] = rect_quad_multipareig(A); t2 = toc;

n2 = size(lambda2,1);
err2 = [];
for j = 1:n2
    mu(1) = lambda2(j,1);
    mu(2) = lambda2(j,2);
    mu(3) = lambda2(j,3);
    M = A{1} + mu(1)*A{2} + mu(2)*A{3} + mu(3)*A{4} + mu(1)^2*A{5} + mu(1)*mu(2)*A{6} + ... 
             + mu(1)*mu(3)*A{7} + mu(2)^2*A{8} + mu(2)*mu(3)*A{9} + mu(3)^2*A{10};
    err2(j,1) = norm(M*X2(:,j));
end

fprintf('k=3, matrices (%d x %d), found %d eigenvalues, time: %7.1e, max residual: %7.1e\n',n+2,n,n2,t2,max(err2));

% Four parameters, k = 4
n = 3;
m = 15; 
A = cell(m,1);
for j = 1:m
    A{j} = rand(n+3,n);
end

tic; [lambda3,X3] = rect_quad_multipareig(A); t3 = toc;

n3 = size(lambda3,1);
err3 = [];
for j = 1:n3
    mu(1) = lambda3(j,1);
    mu(2) = lambda3(j,2);
    mu(3) = lambda3(j,3);
    mu(4) = lambda3(j,4);
    M = A{1} + mu(1)*A{2} + mu(2)*A{3} + mu(3)*A{4} + mu(4)*A{5} + mu(1)^2*A{6} + mu(1)*mu(2)*A{7} + ... 
             + mu(1)*mu(3)*A{8} + mu(1)*mu(4)*A{9} + mu(2)^2*A{10} + mu(2)*mu(3)*A{11} + ...
             + mu(2)*mu(4)*A{12} + mu(3)^2*A{13} + mu(3)*mu(4)*A{14} + mu(4)^2*A{15};
    err3(j,1) = norm(M*X3(:,j));
end

fprintf('k=4, matrices (%d x %d), found %d eigenvalues, time: %7.1e, max residual: %7.1e\n',n+3,n,n3,t3,max(err3));
