function [lambda,mu,x,y,res,flag,iter] = twopareig_refine(A1,B1,C1,A2,B2,C2,lambda,mu,x0,y0,maxiter,tol)

%TWOPAREIG_REFINE  Newton refinement for a two-parameter eigenvalue problem
%
% [lambda,mu,x,y,res,flag,iter] = TWOPAREIG_REFINE(A1,B1,C1,A2,B2,C2,lambda0,mu0,x0,y0,maxiter,tol)
% applies Newton iteration to a two-parameter eigenvalue problem
%
% A1 x = lambda B1 x + mu C1 x
% A2 y = lambda B2 y + mu C2 y
%
% Input
%   - A1,B1,C1,A2,B2,C2 : matrices of the two-parameter eigenvalue problem
%   - lambda, mu : approximation for the eigenvalue
%   - x0, y0 : approximation for the eigenvector
%   - maxiter : maximum number of iterations (default is 1)
%   - tol : tolerance for the residuals (default is 1e2*eps)
%
% Output:
%   - lambda, mu : eigenvalue
%   - x,y : eigenvector
%   - res : norms of the residual
%   - flag : convergence flag (success 2, improvement 1 or fail 0)
%   - iter : number of iterations
%
% See also: MULTIPAREIG_REFINE

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 22.10.2018
% BP 23.11.2023: use relative tolerance

% This file is kept for backward compatibility, use multipareig_refine instead

% Validate number of input parameters
narginchk(10,12);

class_t = superiorfloat(A1,B1,C1,A2,B2,C2,lambda,mu,x0,y0);

% Make sure all inputs are of the same numeric type.
if ~isa(A1,class_t), A1 = numeric_t(A1,class_t); end
if ~isa(B1,class_t), B1 = numeric_t(B1,class_t); end
if ~isa(C1,class_t), C1 = numeric_t(C1,class_t); end
if ~isa(A2,class_t), A2 = numeric_t(A2,class_t); end
if ~isa(B2,class_t), B2 = numeric_t(B2,class_t); end
if ~isa(C2,class_t), C2 = numeric_t(C2,class_t); end
if ~isa(lambda,class_t), lambda = numeric_t(lambda,class_t); end
if ~isa(mu,class_t), mu = numeric_t(mu,class_t); end
if ~isa(x0,class_t), x0 = numeric_t(x0,class_t); end
if ~isa(y0,class_t), y0 = numeric_t(y0,class_t); end

if nargin<11 || isempty(maxiter), maxiter = 1; end
if nargin<12 || isempty(tol), tol = numeric_t('1e2*eps',class_t); end

maxnorm1 = norm(A1,'fro') + abs(lambda)*norm(B1,'fro') + abs(mu)*norm(C1,'fro');
maxnorm2 = norm(A2,'fro') + abs(lambda)*norm(B2,'fro') + abs(mu)*norm(C2,'fro');
tol = tol*max(maxnorm1,maxnorm2);

iter = 1;

x = x0/norm(x0);
y = y0/norm(y0);

n1 = length(x0);
n2 = length(y0);
N = n1+n2;

A1x = A1*x; B1x = B1*x; C1x = C1*x;
A2y = A2*y; B2y = B2*y; C2y = C2*y;

W1 = A1 - lambda*B1 - mu*C1;
W2 = A2 - lambda*B2 - mu*C2;
vec1 = A1x - lambda*B1x - mu*C1x;
vec2 = A2y - lambda*B2y - mu*C2y;
res(1) = norm(vec1);
res(2) = norm(vec2);
residual(iter) = norm(res);

while iter<= maxiter && residual(iter) >= tol
    
    J11 = [W1 zeros(n1,n2); zeros(n2,n1) W2];
    J12 = [B1x C1x; B2y C2y];
    J21 = [x0' zeros(1,n2); zeros(1,n1) y0'];
    J22 = zeros(2, class_t);
    
    % build the Jacobian matrix and the right hand side for the Newton step
    J = [J11 -J12; J21 J22];
    RH = [vec1; vec2; x0'*x-1; y0'*y-1];
    
    warning off
    delta = -J\RH;
    warning on
    x = x + delta(1:n1,1);
    y = y + delta(n1+1:n1+n2,1);
    lambda = lambda + delta(end-1,1);
    mu = mu + delta(end,1);
    
    A1x = A1*x; B1x = B1*x; C1x = C1*x;
    A2y = A2*y; B2y = B2*y; C2y = C2*y;

    W1 = A1 - lambda*B1 - mu*C1;
    W2 = A2 - lambda*B2 - mu*C2;
    vec1 = A1x - lambda*B1x - mu*C1x;
    vec2 = A2y - lambda*B2y - mu*C2y;
    res(1) = norm(vec1);
    res(2) = norm(vec2);
    
    iter = iter + 1;
    residual(iter) = norm(res);
end

x = x/norm(x);
y = y/norm(y);

if residual(iter)<tol 
    flag = 2;           
elseif residual(iter)<residual(1)
    flag = 1;  
else
    flag = 0; 
end
