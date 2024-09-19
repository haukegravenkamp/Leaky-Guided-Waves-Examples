function [lambda,mu,X] = rect_quad_twopareig(A,B,C,D,E,F,opts)

%RECT_QUAD_TWOPAREIG  Solve a rectangular quadratic two-parameter eigenvalue problem
%
% [lambda,mu,X] = RECT_QUAD_TWOPAREIG(A,B,C) returns eigenvalues and 
% eigenvectors of the rectangular quadratic two-parameter eigenvalue problem
%
% (A + lambda*B + mu*C + lambda^2*D+ lambda*mu*E + mu^2*F)x = 0
% 
% where A,B,C,D,E,F are (n+1) x n matrices
% 
% Input:
%   - A,B,C,D,E,F : matrices of size (n+1) x n
%   - opts : options
%
% Options in opts:
%   - method : 'compress' (default) or 'mep' or 'linearize'
%   - delta (sqrt(eps)): treshold for eigenvalues of the rectangular problem
%   - all options for related method quad_twopareig for 'mep' or
%     joint_delta_eig for ('compress','linear')
% Output:
%   - lambda,mu : eigenvalues
%   - X : matrix of size n x m with right eigenvectors
%
% This file is kept for backward compatibility, use rect_quad_multipareig instead

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% P. Holoborodko, Advanpix LLC.
% FreeBSD License, see LICENSE.txt

% Bor Plestenjak 2022

% See: M.E.Hochstenbach, T.Kosir, B.Plestenjak: Numerical methods for rectangular 
% multiparameter eigenvalue problems, with applications to finding optimal 
% ARMA and LTI models. Numer Linear Algebra Appl. 2023; e2540

% Compression method is based on: F. F. Alsubaie: H2 Optimal Model 
% Reduction for Linear Dynamic Systems and the Solution of Multiparameter 
% Matrix Pencil Problems, PhD, Imperial College London, 2019.


narginchk(6,7);

% Analyse user supplied options, if any
if nargin < 7, opts = []; end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(A,B,C,D,E,F);
end

if isfield(opts,'method'),   method  = opts.method;   else,   method = 'mep';          end
if isfield(opts,'delta'),    delta   = opts.delta;    else,   delta = sqrt(eps(class_t));   end

if ~isa(A,class_t), A = numeric_t(A,class_t); end
if ~isa(B,class_t), B = numeric_t(B,class_t); end
if ~isa(C,class_t), C = numeric_t(C,class_t); end
if ~isa(D,class_t), D = numeric_t(D,class_t); end
if ~isa(E,class_t), E = numeric_t(E,class_t); end
if ~isa(F,class_t), F = numeric_t(F,class_t); end

[m,n] = size(A);
if m ~= n+1
   error('Matrices must be rectangular of size (n+1) x n')
end
[m1,n1] = size(B);
[m2,n2] = size(C);
[m3,n3] = size(D);
[m4,n4] = size(E);
[m5,n5] = size(F);
if any([m1 n1] ~= [m n]) || any([m2 n2] ~= [m n]) || any([m3 n3] ~= [m n]) || any([m4 n4] ~= [m n]) || any([m5 n5] ~= [m n])
   error('Matrices must be rectangular of size (n+1) x n')
end

if strcmp(method,'compress')
    if nargout>2
        [Lambda,X] = rect_quad_multipareig({A,B,C,D,E,F},opts);
    else
        Lambda = rect_quad_multipareig({A,B,C,D,E,F},opts);
    end
    lambda = Lambda(:,1);
    mu = Lambda(:,2);
elseif strcmp(method,'mep')
    % random projection matrices
    P1 = orth(randn(m,n,class_t)).';
    P2 = orth(randn(m,n,class_t)).';
    % we project original problem into a two-parameter eigenvalue problem
    A1 = P1*A; B1 = P1*B; C1 = P1*C; D1 = P1*D; E1 = P1*E; F1 = P1*F;
    A2 = P2*A; B2 = P2*B; C2 = P2*C; D2 = P2*D; E2 = P2*E; F2 = P2*F;

    nA = norm(A,'fro');
    nB = norm(B,'fro');
    nC = norm(C,'fro');
    nD = norm(D,'fro');
    nE = norm(E,'fro');
    nF = norm(F,'fro');

    % we solve the two-parameter eigenvalue problem that has extra solutions
    % and then return only solutions that solve the rectangular problem
    [lambdaC,muC] = quad_twopareig(A1,B1,C1,D1,E1,F1,A2,B2,C2,D2,E2,F2,opts);
    lambda = [];
    mu = [];
    X = [];
    for j=1:length(lambdaC)
      mat = A + lambdaC(j)*B + muC(j)*C + lambdaC(j)^2*D + lambdaC(j)*muC(j)*E + muC(j)^2*F;
      ocena = nA + abs(lambdaC(j))*nB + abs(muC(j))*nC + abs(lambdaC(j))^2*nD + abs(lambdaC(j))*abs(muC(j))*nE + abs(muC(j))^2*nF;
      xr = min_sing_vec(mat,0);
      res = norm(mat*xr);
      % we select eigenvalues that give small residual for the original problem
      if norm(res)<delta*ocena
          lambda = [lambda; lambdaC(j)];
          mu = [mu; muC(j)];
          X = [X xr];
      end
    end
elseif strcmp(method,'linearize') 
    % linearization to linear rectangular MEP
    ZX = zeros(n+1,n,class_t);
    Z = zeros(n,class_t);
    I = eye(n,class_t);
    RA = sparse([A  B  C; Z  -I  Z; Z  Z  -I]);
    RB = sparse([ZX D  E; I   Z  Z; Z  Z   Z]);
    RC = sparse([ZX ZX F; Z   Z  Z; I  Z   Z]);
   
    [RDelta0,RDelta1,RDelta2] = twopar_delta(RA,-RB,-RC,RA,-RB,-RC);
    % standard compression for linear RMEP
    PTR = right_compression_matrix(3*n,2);
    PTL = left_compression_matrix(3*n+1,2);
    DW0 = PTL*RDelta0*PTR;
    DW1 = PTL*RDelta1*PTR;
    DW2 = PTL*RDelta2*PTR;
    
    opts.singular = 1;
    lamu = joint_delta_eig({DW0,DW1,DW2},opts);
    lambda = lamu(:,1);
    mu = lamu(:,2);
    m = length(lambda);
    
    if nargout>2
        X = zeros(n,m,class_t);
        for j = 1:m
           M = A + lambda(j)*B + mu(j)*C + lambda(j)^2*D + lambda(j)*mu(j)*E + mu(j)^2*F;
           X(:,j) = min_sing_vec(M,0);
        end
    end
else
    error('Unknown method specified in options. Use "mep", "linearize" or "compress"')
end
