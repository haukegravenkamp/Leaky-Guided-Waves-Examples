function [lambda,mu,eta,X1,X2,X3,Y1,Y2,Y3,report] = threepareig(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,opts)

%THREEPAREIG   Solve a three-parameter eigenvalue problem
%
% [lambda,mu,eta,X1,X2,X3,Y1,Y2,Y3] = THREEPAREIG(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,opts) 
% finds eigenvalues and eigenvectors of a three-parameter eigenvalue problem
%
% A1 x1 = lambda B1 x1 + mu C1 x1 + eta D1 x1 
% A2 x2 = lambda B2 x2 + mu C2 x2 + eta D2 x2 
% A3 x3 = lambda B3 x3 + mu C3 x3 + eta D3 x3 
%
% Input:
%   - A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3 : matrices
%   - opts : options (see MULTIPAREIG)
%
% Output:
%   - lambda , mu, eta: eigenvalues (eigenvalue is (lambda(j),mu(j),eta(j))
%   - X1, X2, X3 : components of decomposable right eigenvectors 
%     (eigenvector is kron(X1(:,j),kron(X2(:,j),X3(:,j))), such that
%       (A1-lambda(j)*B1-mu(j)*C1-eta(j)*D1)*X1(:,j)=0
%       (A2-lambda(j)*B2-mu(j)*C2-eta(j)*D2)*X2(:,j)=0
%       (A3-lambda(j)*B3-mu(j)*C3-eta(j)*D3)*X3(:,j)=0
%   - Y1, Y2, Y3 : components of decomposable left eigenvectors 
%     (eigenvector is kron(Y1(:,j),kron(Y2(:,j),Y3(:,j))), such that
%       (A1-lambda(j)*B1-mu(j)*C1-eta(j)*D1)'*Y1(:,j)=0
%       (A2-lambda(j)*B2-mu(j)*C2-eta(j)*D2)'*Y2(:,j)=0
%       (A3-lambda(j)*B3-mu(j)*C3-eta(j)*D3)'*Y3(:,j)=0
%
% Operator determinants Delta0, Delta1, Delta2, Delta3 are used, where
% Delta0 =   | B1 C1 D1; B2 C2 D2; B3 C3 D3 |
% Delta1 = - | A1 C1 D1; A2 C2 D2; A3 C3 D3 |
% Delta2 =   | B1 A1 D1; B2 A2 D2; B3 A3 D3 |
% Delta3 = - | B1 C1 A1; B2 C2 A2; B3 C3 A3 |
%
% See also: TWOPAREIG, MULTIPAREIG, THREEPAREIGS, THREEPAREIGS_JD,
% THREEPAREIGS_SI

% Reference: M. E. Hochstenbach, T. Kosir, B. Plestenjak: A Jacobi-Davidson 
% type method for the two-parameter eigenvalue problem, SIAM J. Matrix Anal. 
% Appl. 26 (2005) 477-497

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% P. Holoborodko, Advanpix LLC.
% FreeBSD License, see LICENSE.txt

% BP 06.09.2015 : support for singular 3EP
% BP 03.11.2016 : small speedup in computation of mu and eta (inspired by Pavel Holoborodko's changes in multipareig)
% PH 22.11.2016 : fixed error when fast=0, added precision-independency.
% BP 26.11.2016 : option fp_type
% PH 26.11.2016 : code simplifications and clean-ups.
% BP 03.04.2022 : computation is moved into multipareig
% Last revision: 03.04.2022

% This file is kept for backward compatibility, use multipareig instead

% Validate number of input parameters.
narginchk(12, 13);

if nargin < 13, opts = []; end
if isfield(opts,'novectors'),  novectors  = opts.novectors;  else, novectors  = 0;  end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3);
end

if nargout<4 
    novectors = 1; 
end

% Default outputs
lambda = numeric_t([],class_t);  %#ok<*NASGU>
eta    = numeric_t([],class_t); 
mu     = numeric_t([],class_t); 
X1     = numeric_t([],class_t);  
X2     = numeric_t([],class_t);  
X3     = numeric_t([],class_t);  
Y1     = numeric_t([],class_t);  
Y2     = numeric_t([],class_t);  
Y3     = numeric_t([],class_t);  
report = numeric_t([],class_t);

% we call multipareig
A = {A1,B1,C1,D1; A2,B2,C2,D2; A3,B3,C3,D3};
opts.novectors = novectors;
[Lambda,X,Y,report] = multipareig(A,opts);

m = size(Lambda,1);
if m>0
    lambda = Lambda(:,1);
    mu     = Lambda(:,2);
    eta    = Lambda(:,3);
    if ~novectors
        X1 = zeros(size(A1,1),m,class_t);
        X2 = zeros(size(A2,1),m,class_t);        
        X3 = zeros(size(A3,1),m,class_t);
        Y1 = zeros(size(A1,1),m,class_t);
        Y2 = zeros(size(A2,1),m,class_t);        
        Y3 = zeros(size(A3,1),m,class_t);
        for j = 1:m
            X1(:,j) = X{j,1};  X2(:,j) = X{j,2};  X3(:,j) = X{j,3};  
            Y1(:,j) = Y{j,1};  Y2(:,j) = Y{j,2};  Y3(:,j) = Y{j,3};  
        end
    end
end

