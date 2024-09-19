function [lambda,mu,eta,X,Y,Z] = threepareigs(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,k,opts)

%THREEPAREIGS   Find a few eigenvalues for a three-parameter eigenvalue problem
%
% [lambda,mu,eta,X,Y,Z] = THREEPAREIGS(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,neig,opts)
% returns neig eigenvalues with the smallest |eta| of the three-parameter
% eigenvalue problem
%
% A1*x = lambda*B1*x + mu*C1*x + eta*D1*x
% A2*y = lambda*B2*y + mu*C2*y + eta*D2*y
% A3*z = lambda*B3*z + mu*C3*z + eta*D3*z
%
% using implicit restarted inverse Arnoldi in Matlab routine eigs
% on the related generalized eigenvalue problems Delta3*w = eta*Delta0*w, 
% where Delta0 and Delta3 are corresponding operator determinants
% Delta0 =   | B1 C1 D1; B2 C2 D2; B3 C3 D3 |
% Delta3 = - | B1 C1 A1; B2 C2 A2; B3 C3 A3 |
%
% Input:
%   - A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3 : matrices
%   - k : number of eigenvalues
%   - opts : options (see below)
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
% Options in opts:
%   - usesparse : set to 0 if all matrices are full. The default (1) uses 
%     sparse representation of Delta matrices, which is better if some 
%     matrices are diagonal or zero
%   - fp_type: numeric type to use ('single', 'double', or 'mp' (needs MCT),
%     use only if you want to change the default - the superior type of input data
%   - refine : number of TRQ refinement steps in the end (2)
%   - all options for the eigs
%
% See also: THREEPAREIG, THREEPAREIGS_JD, THREEPAREIGS_SI, TWOPAREIGS

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 26.08.2015 : new option usesparse
% BP 03.09.2015 : refine with trqi_3p
% BP 13.09.2016 : efficient computation of Rayleigh quotients
% BP 01.02.2018 : modified to be precision-independent 
% BP 26.02.2018 : more efficient method - we need only two Delta matrices

% Last revision: 01.02.2018

narginchk(12,14);

% Analyse user supplied options, if any.
if nargin < 13, k = 6; end   % if not specified, 6 eigenvalues are returned
if nargin < 14, opts = []; end
if isfield(opts,'usesparse'),  usesparse = opts.usesparse;  else usesparse = 1;  end
if isfield(opts,'refine'),     refine = opts.refine;  else refine = 2;           end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3);
end

% Make sure all inputs are of the same numeric type.
if ~isa(A1,class_t), A1 = numeric_t(A1,class_t); end;
if ~isa(B1,class_t), B1 = numeric_t(B1,class_t); end;
if ~isa(C1,class_t), C1 = numeric_t(C1,class_t); end;
if ~isa(D1,class_t), D1 = numeric_t(D1,class_t); end;
if ~isa(A2,class_t), A2 = numeric_t(A2,class_t); end;
if ~isa(B2,class_t), B2 = numeric_t(B2,class_t); end;
if ~isa(C2,class_t), C2 = numeric_t(C2,class_t); end;
if ~isa(D2,class_t), D2 = numeric_t(D2,class_t); end;
if ~isa(A3,class_t), A3 = numeric_t(A3,class_t); end;
if ~isa(B3,class_t), B3 = numeric_t(B3,class_t); end;
if ~isa(C3,class_t), C3 = numeric_t(C3,class_t); end;
if ~isa(D3,class_t), D3 = numeric_t(D3,class_t); end;

n1 = size(A1,1);
n2 = size(A2,1);
n3 = size(A3,1);

% if k>=n1*n2*n3 we compute all eigenvalues 
if k>=n1*n2*n3-1
    [lambda,mu,eta,X,Y,Z] = threepareig(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3);
    return
end

if usesparse 
    A1s = sparse(A1); B1s = sparse(B1); C1s = sparse(C1); D1s = sparse(D1);
    A2s = sparse(A2); B2s = sparse(B2); C2s = sparse(C2); D2s = sparse(D2);
    A3s = sparse(A3); B3s = sparse(B3); C3s = sparse(C3); D3s = sparse(D3);
    Delta0 = kron(B1s,kron(C2s,D3s)-kron(D2s,C3s)) - kron(C1s,kron(B2s,D3s)-kron(D2s,B3s)) + kron(D1s,kron(B2s,C3s) - kron(C2s,B3s));
    Delta3 = kron(B1s,kron(C2s,A3s)-kron(A2s,C3s)) - kron(C1s,kron(B2s,A3s)-kron(A2s,B3s)) + kron(A1s,kron(B2s,C3s) - kron(C2s,B3s));
else
    Delta0 = kron(B1,kron(C2,D3)-kron(D2,C3)) - kron(C1,kron(B2,D3)-kron(D2,B3)) + kron(D1,kron(B2,C3) - kron(C2,B3));
    Delta3 = kron(B1,kron(C2,A3)-kron(A2,C3)) - kron(C1,kron(B2,A3)-kron(A2,B3)) + kron(A1,kron(B2,C3) - kron(C2,B3));
end

opts.disp = 0; % turn off diagnostic information level in eigs in earlier versions of Matlab
[W,D] = eigs(Delta3,Delta0,k,'SM',opts);
eta = diag(D);

% lambda and mu are computed from Rayleigh quotients
for i=1:k
    if usesparse
        [y0,y1,y2,y3] = threepar_delta_mv(W(:,i),A1s,B1s,C1s,D1s,A2s,B2s,C2s,D2s,A3s,B3s,C3s,D3s);
    else
        [y0,y1,y2,y3] = threepar_delta_mv(W(:,i),A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3);
    end
    imen = (W(:,i)'*y0);
    mu(i,1) = W(:,i)'*y2/imen;
    lambda(i,1) = W(:,i)'*y1/imen;
end

if nargout > 3
    % extraction of eigenvectors (individually using inverse iteration or SVD)
    for j=1:k
        [xr,xl] = min_sing_vec(A1-lambda(j)*B1-mu(j)*C1-eta(j)*D1);
        X(:,j) = xr;
        [xr,xl] = min_sing_vec(A2-lambda(j)*B2-mu(j)*C2-eta(j)*D2);
        Y(:,j) = xr;
        [xr,xl] = min_sing_vec(A3-lambda(j)*B3-mu(j)*C3-eta(j)*D3);
        Z(:,j) = xr;
        if refine>0
            [refl,refu,refe,refx,refy,refz] = trqi_3p(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,X(:,j),Y(:,j),Z(:,j),refine,eps);
            X(:,j) = refx;  Y(:,j) = refy;  Z(:,j) = refz;
            lambda(j) = refl;  mu(j) = refu;  eta(j) = refe;
        end
    end   
end

