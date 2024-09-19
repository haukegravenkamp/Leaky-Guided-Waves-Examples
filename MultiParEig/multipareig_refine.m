function [lambda,X,res,flag,iter] = multipareig_refine(A,lambda,X0,maxiter,tol)

%MULTIPAREIG_REFINE  Newton refinement for a multiparameter eigenvalue problem
%
% [lambda,X,res,flag,iter] = MULTIPAREIG_REFINE(A,lambda,X,maxiter,tol)
% applies Newton iteration to a multiparameter eigenvalue problem
%
% A{1,1} x1 = lambda(1) A{1,2} x1 +  ... + lambda(k) A{1,k+1} x1 
% A{2,1} x2 = lambda(1) A{2,2} x2 +  ... + lambda(k) A{2,k+1} x2 
% ...
% A{k,1} xk = lambda(1) A{k,2} xk +  ... + lambda(k) A{k,k+1} xk 
%
% Input:
%   - A : cell array of size k x (k+1) of matrices A{i,j}, all matrices in
%         the same row have to be square matrices of the same size
%   - lambda : eigenvalue, row of size k
%   - X : eigenvector, cell array of size 1 x k with right eigenvectors
%   - maxiter : maximum number of iterations (default is 1)
%   - tol : tolerance for the residuals (default is 1e2*eps)
%
% Output:
%   - lambda : eigenvalue, row of size k
%   - X : eigenvector as cell arrray
%   - res : norms of the residual
%   - flag : convergence flag (success 2, improvement 1, fail 0)
%   - iter : number of iterations
%
% See also: TWOPAREIG_REFINE

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 23.11.2023: use relative tolerance

% Validate number of input parameters
narginchk(3,5);

class_t = superiorfloat(A{:});

% Make sure all inputs are of the same numeric type.
for k=1:numel(A)
    if ~isa(A{k}, class_t)
         A{k} = numeric_t(A{k},class_t);
    end
end
if ~isa(lambda,class_t), lambda = numeric_t(lambda,class_t); end
for k=1:numel(X0)
    if ~isa(X0{k}, class_t)
         X0{k} = numeric_t(X0{k},class_t);
    end
end

maxnorm = zeros(k,1);
for i = 1:k
    maxnorm(i) = norm(A{i,1},'fro');
    for j = 1:k
        maxnorm(i) = maxnorm(i) + abs(lambda(j))*norm(A{i,j+1},'fro');
    end
end

if nargin<4 || isempty(maxiter), maxiter = 1; end
if nargin<5 || isempty(tol), tol = numeric_t('1e2*eps',class_t); end
tol = tol*norm(maxnorm,"inf");

k = size(A,1);
N = 0;
left = 1;
right = 0;
for j = 1:k
    X{j} = X0{j}/norm(X0{j});
    n(j) = size(A{j,1},1);
    nl(j) = left;
    left = left + n(j);
    right = right + n(j);
    nr(j) = right;
end
N = sum(n);

AX = cell(k,k+1);
J11 = zeros(N,class_t);
J12 = zeros(N,k,class_t);
J21 = zeros(k,N,class_t);
J22 = zeros(k,class_t);
RH = zeros(N+k,1,class_t);
vec = cell(k,1);

iter = 1;
for i = 1:k
    for j = 1:k+1
        AX{i,j} = A{i,j}*X{i};
    end
end

for i = 1:k
    Wi = A{i,1};
    vec{i} = AX{i,1};
    for j = 1:k
        Wi = Wi - lambda(j)*A{i,j+1};
        J12(nl(i):nr(i),j) = AX{i,j+1};
        vec{i} = vec{i} - lambda(j)*AX{i,j+1};
    end
    J11(nl(i):nr(i),nl(i):nr(i)) = Wi;
    J21(i,nl(i):nr(i)) = X0{i}';
    RH(nl(i):nr(i)) = vec{i};
    RH(N+i) = X0{i}'*X{i}-1;
    res(i) = norm(vec{i});
end
residual(iter) = norm(res);

while iter<= maxiter && residual(iter)>=tol

    % build the Jacobian matrix and the right hand side for the Newton step
    J = [J11 -J12; J21 J22];
    
    warning off
    delta = -J\RH;
    warning on
    for i = 1:k
        X{i} = X{i} + delta(nl(i):nr(i),1);
    end
    lambda = lambda + delta(N+1:end,1).';

    % build matrices for the next step and compute residual
    for i = 1:k
        for j = 1:k+1
            AX{i,j} = A{i,j}*X{i};
        end
    end
    
    for i = 1:k
        Wi = A{i,1};
        vec{i} = AX{i,1};
        for j = 1:k
            Wi = Wi - lambda(j)*A{i,j+1};
            J12(nl(i):nr(i),j) = AX{i,j+1};
            vec{i} = vec{i} - lambda(j)*AX{i,j+1};
        end
        J11(nl(i):nr(i),nl(i):nr(i)) = Wi;
        J21(i,nl(i):nr(i)) = X0{i}';
        RH(nl(i):nr(i)) = vec{i};
        RH(N+i) = X0{i}'*X{i}-1;
        res(i) = norm(vec{i});
    end
    iter = iter + 1;
    residual(iter) = norm(res);
end

for i=1:k
    X{i} = X{i}/norm(X{i});
end

if norm(res)<tol 
    flag = 2;           
elseif residual(iter)<residual(1)
    flag = 1;  
else
    flag = 0; 
end

