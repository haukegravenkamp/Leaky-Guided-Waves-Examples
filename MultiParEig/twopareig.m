function [lambda,mu,XR,YR,XL,YL,report] = twopareig(A1,B1,C1,A2,B2,C2,opts)

%TWOPAREIG   Solve a two-parameter eigenvalue problem
%
% [lambda,mu,XR,YR,XL,YL] = TWOPAREIG(A1,B1,C1,A2,B2,C2,opts) returns
% eigenvalues and eigenvectors of the two-parameter eigenvalue problem
%
% A1 x = lambda B1 x + mu C1 x 
% A2 y = lambda B2 y + mu C2 y
%
% Input:
%   - A1, B1, C1, A2, B2, C2: matrices
%   - opts: options (see below)
%
% Output: 
%   - lambda, mu: eigenvalue parts (eigenvalues are (lambda(j),mu(j))
%   - XR, YR: components of decomposable right eigenvectors
%     (eigenvector is kron(XR(:,j),YR(:,j)) such that 
%       (A1-lambda(j)*B1-mu(j)*C1)*XR(:,j)=0
%       (A2-lambda(j)*B2-mu(j)*C2)*YR(:,j)=0
%   - XL, YL: components of decomposable left eigenvectors 
%     (eigenvector is kron(XL(:,j),YL(:,j)) such that 
%       (A1-lambda(j)*B1-mu(j)*C1)'*XL(:,j)=0
%       (A2-lambda(j)*B2-mu(j)*C2)'*YL(:,j)=0
%   - report: details of compression method in case of a singular problem 
%       rows are [mode m n r s(1) s(r) s(r+1) choice)], 
%       where mode is 1: CR step 1, 2: CR step 2, 3: RC step 1, 4. RC step 2
%       m,n: size of Delta matrices, r: rank
%       s(1), s(r), s(r+1): singular values
%       choice: how was rank determined (see NUMRANK for details)
%
% Operator determinants Delta0, Delta1, and Delta2 are used, where
% Delta0 = kron(C2, B1) - kron(B2, C1)
% Delta1 = kron(C2, A1) - kron(A2, C1)
% Delta2 = kron(A2, B1) - kron(B2, A1)
%
% Options in opts:
%   - novectors (0): set to 1 when report is needed and XR,YR,XL,YL are not
%   - fp_type: numeric type to use ('single', 'double', or 'mp' (needs MCT),
%     use only if you want to change the default - the superior type of input data
%   - inviter (1): use inverse iteration for eigenvectors or slow svd (0)
%   - refine (1): Newton refinement steps to improve the accuracy of simple
%     eigenvalues of a regular MEP - requires computation of eigenvectors 
%   - refineeps (eps): relative tolerance for Newton refinement
% Options in opts related to JOINT_DELTA_EIG:
%   - singular (0): set to 1 for a singular problem, i.e., det(Delta0)=0
%   - rand_orth (1): take DeltaH as a random combination of Delta1 and Delta2,
%     set to 0 to take DeltaH=Delta1
%   - solver ('eig'): numerical approach to find eigenvalues:
%       - 'eig': compute eigenvectors of inv(Delta0)*DeltaH and then
%          eigenvalues of inv(Delta0)*Deltaj from Rayleigh quotients
%       - 'geneig': compute eigenvectors of (DeltaH, Delta0) and then
%          eigenvalues of (Deltaj,Delta0) from Rayleigh quotients
%       - 'schur': compute Schur form of inv(Delta0)*DeltaH with optional
%          clustering and reordering and use it to (block) triangularize
%          inv(Delta0)*Deltaj for j=1,2
%       - 'genschur': compute generalized Schur form of (DeltaH,Delta0) 
%          with optional clustering and reordering and use it to (block) 
%          triangularize (Deltaj,Delta0) for j=1,2
%   - twosidedRQ (0): set to 1 to use two-sided Rayleigh quotients instead
%     of one-sided (only for solvers 'eig' and 'geneig') 
%   - epscluster (1e-6): relative distance between eigenvalues in a cluster,
%     set to 0 for no clustering (only for solvers 'schur' and 'genschur')
%   - rrqr (0): set to 1 to use rank revealing QR instead of SVD in the
%     staircase-type method (only for singular problems)
%   - rng_seed (0): if different from zero, rng(rng_seed) is used at start
%   - epsshape (1e-6): displays warning if shape error of simultaneous 
%     block triangularization is larger that epsshape, set to 0 for no testing

% See also: MULTIPAREIG, TWOPAREIGS_JD, TWOPAREIGS_SI, THREEPAREIG

% Reference: M. E. Hochstenbach, T. Kosir, B. Plestenjak: A Jacobi-Davidson 
% type method for the two-parameter eigenvalue problem, SIAM J. Matrix Anal. 
% Appl. 26 (2005) 477-497

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% P. Holoborodko, Advanpix LLC.
% FreeBSD License, see LICENSE.txt

% BP 06.09.2015 : use extract_regular_part_np
% BP 24.08.2016 : double clustering in QZ (lambda first and then mu in blocks)
% BP 03.11.2016 : small speedup in computation od mu (inspired by Pavel Holoborodko's changes in multipareig)
% PH 22.11.2016 : modified to be precision-independent, i.e. be able to work with 
%                 any numeric type (whether it is built-in 'double'/'single' or custom like 'mp').
% BP 26.11.2016 : option fp_type
% PH 26.11.2016 : code simplifications, small fixes and clean-ups.
% BP 22.10.2018 : Newton refinement
% BP 03.04.2022 : random combination for joint triangularization
% BP 11.11.2023 : computation is moved into multipareig
% BP 16.01.2024 : cleaning of options

% This file is kept for backward compatibility, use multipareig instead

narginchk(6, 7);

% Analyse user supplied options, if any
if nargin < 7, opts = []; end
if isfield(opts,'novectors'),  novectors  = opts.novectors;  else, novectors  = 0;  end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(A1,B1,C1,A2,B2,C2);
end

% Make sure all inputs are of the same numeric type.
if ~isa(A1,class_t), A1 = numeric_t(A1,class_t); end
if ~isa(B1,class_t), B1 = numeric_t(B1,class_t); end
if ~isa(C1,class_t), C1 = numeric_t(C1,class_t); end
if ~isa(A2,class_t), A2 = numeric_t(A2,class_t); end
if ~isa(B2,class_t), B2 = numeric_t(B2,class_t); end
if ~isa(C2,class_t), C2 = numeric_t(C2,class_t); end

if nargout<3 
    novectors = 1; 
end

% Default outputs
lambda = numeric_t([],class_t); 
mu     = numeric_t([],class_t); 
XR     = numeric_t([],class_t);  
YR     = numeric_t([],class_t);  
XL     = numeric_t([],class_t);  
YL     = numeric_t([],class_t); 
report = numeric_t([],class_t);

% we call multipareig
A = {A1,B1,C1; A2,B2,C2};
opts.novectors = novectors;
[Lambda,X,Y,report] = multipareig(A,opts);

m = size(Lambda,1);
if m>0
    lambda = Lambda(:,1);
    mu     = Lambda(:,2);
    if ~novectors
        XR = zeros(size(A1,1),m,class_t);
        YR = zeros(size(A2,1),m,class_t);        
        XL = zeros(size(A1,1),m,class_t);
        YL = zeros(size(A2,1),m,class_t);        
        for j = 1:m
            XR(:,j) = X{j,1};  YR(:,j) = X{j,2};  
            XL(:,j) = Y{j,1};  YL(:,j) = Y{j,2};
        end
    end
end