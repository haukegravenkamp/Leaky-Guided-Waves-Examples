function [lambda,X,Y,report] = multipareig(A,opts)

%MULTIPAREIG   Solve a multiparameter eigenvalue problem
%
% [lambda,X,Y] = MULTIPAREIG(A,opts) returns eigenvalues and eigenvectors 
% of the multiparameter eigenvalue problem
%
% A{1,1} x1 = lambda(1) A{1,2} x1 +  ... + lambda(k) A{1,k+1} x1 
% A{2,1} x2 = lambda(1) A{2,2} x2 +  ... + lambda(k) A{2,k+1} x2 
% ...
% A{k,1} xk = lambda(1) A{k,2} xk +  ... + lambda(k) A{k,k+1} xk 
%
% Input:
%   - A : cell array of size k x (k+1) of matrices A{i,j}, all matrices in
%         the same row have to be square matrices of the same size
%   - opts : options
%
% Options in opts:
%   - novectors (0): set to 1 when report is needed and X and Y are not
%   - fp_type: numeric type to use ('single', 'double', or 'mp' (needs MCT),
%     use only if you want to change the default - the superior type of input data
%   - inviter (1): use inverse iteration for eigenvectors or slow svd (0)
%   - refine (1): Newton refinement steps to improve the accuracy of simple
%     eigenvalues of a regular MEP - requires computation of eigenvectors 
%   - refineeps (eps): relative tolerance for Newton refinement
% Options in opts related to JOINT_DELTA_EIG:
%   - singular (0): set to 1 for a singular problem, i.e., det(Delta{1})=0
%   - rand_orth (1): take DeltaH as a random combination of Delta{2} to Delta{k+1}, 
%     set to 0 to take DeltaH=Delta{2}
%   - solver ('eig'): numerical approach to find eigenvalues:
%       - 'eig': compute eigenvectors of inv(Delta{1})*DeltaH and then
%          eigenvalues of inv(Delta{1})*Delta{j} from Rayleigh quotients
%       - 'geneig': compute eigenvectors of (DeltaH, Delta{1}) and then
%          eigenvalues of (Delta{j},Delta{1}) from Rayleigh quotients
%       - 'schur': compute Schur form of inv(Delta{1})*DeltaH with optional
%          clustering and reordering and use it to (block) triangularize
%          inv(Delta{1})*Delta{j} for j=2:k+1
%       - 'genschur': compute generalized Schur form of (DeltaH,Delta{1}) 
%          with optional clustering and reordering and use it to (block) 
%          triangularize (Delta{j},Delta{1}) for j=2:k+1
%   - twosidedRQ (0): set to 1 to use two-sided Rayleigh quotients instead
%     of one-sided (only for solvers 'eig' and 'geneig') 
%   - epscluster (1e-6): relative distance between eigenvalues in a cluster,
%     set to 0 for no clustering (only for solvers 'schur' and 'genschur')
%   - rrqr (0): set to 1 to use rank revealing QR instead of SVD in the
%     staircase-type method (only for singular problems)
%   - rng_seed (0): if different from zero, rng(rng_seed) is used at start
%   - epsshape (1e-6): displays warning if shape error of simultaneous 
%     block triangularization is larger that epsshape, set to 0 for no testing
%
% Output:
%   - lambda : matrix of size m x k, each row is an eigenvalue
%   - X : cell array of size m x k with right eigenvectors
%   - Y : cell array of size m x k with left eigenvectors
%   - report: details of compression method in case of a singular problem 
%       rows are [mode m n r s(1) s(r) s(r+1) choice)], 
%       where mode is 1: CR step 1, 2: CR step 2, 3: RC step 1, 4. RC step 2
%       m,n: size of Delta matrices, r: rank
%       s(1), s(r), s(r+1): singular values
%       choice: how was rank determined (see NUMRANK for details)
% 
% See also: TWOPAREIG, THREEPAREIG, JOINT_DELTA_EIG, MULTIPAR_DELTA.

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% P. Holoborodko, Advanpix LLC.
% FreeBSD License, see LICENSE.txt

% BP 05.09.2015 : support for singular MEP
% PH 01.11.2016 : e-values computation speed-up
% PH 22.11.2016 : precision-independent version 
% BP 26.11.2016 : option fp_type
% PH 26.11.2016 : code simplifications and clean-ups.
% BP 02.04.2022 : random combination for joint triangularization
% BP 03.07.2022 : correct sign of delta matrices
% BP 06.12.2022 : eigs of commuting matrices are now in auxiliary file
% BP 17.12.2023 : eigenvectors are not computed if not needed
% BP 16.01.2024 : cleaning of options

% Validate number of input parameters
narginchk(1, 2);

% Analyse user supplied options, if any.
if nargin < 2, opts = []; end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(A{:});
end

if isfield(opts,'inviter'),    inviter    = opts.inviter;    else, inviter    = 1;                         end
if isfield(opts,'novectors'),  novectors  = opts.novectors;  else, novectors  = 0;                         end
if isfield(opts,'refine'),     refine     = opts.refine;     else, refine     = 1;                         end
if isfield(opts,'refineeps'),  refineeps  = opts.refineeps;  else, refineeps  = eps(class_t);              end

% Make sure all inputs are of the same numeric type.
for k=1:numel(A)
    if ~isa(A{k}, class_t)
         A{k} = numeric_t(A{k},class_t);
    end
end

if nargout<2 
    novectors = 1; 
end

% Default outputs
lambda = numeric_t([],class_t); 
X      = numeric_t([],class_t);  
Y      = numeric_t([],class_t);  
report = numeric_t([],class_t);

k = length(A); % number of parameters + 1

% computation of Delta matrices
Delta = multipar_delta(A);

[lambda,report] = joint_delta_eig(Delta,opts);

% extraction of eigenvectors (individually using inverse iteration)
if ((~novectors) && (nargout > 1)) || (refine>0)
    m = size(lambda,1);
    X = cell(m,k-1);
    Y = cell(m,k-1);
    for j = 1:m
        for r = 1:k-1
            M = A{r,1};
            for p = 1:k-1
                M = M - lambda(j,p)*A{r,p+1};
            end
            [zr,zl] = min_sing_vec(M,inviter);
            X{j,r} = zr;
            Y{j,r} = zl;
        end
        if refine>0
            W = cell(1,k-1);
            for p = 1:k-1
                W{p} = X{j,p};
            end
            [newlambda,W1,~,flag] = multipareig_refine(A,lambda(j,:),W,refine,refineeps);
            if flag
                lambda(j,:) = newlambda;
                for p = 1:k-1
                    X{j,p} = W1{p};
                end
            end
        end
    end   
end

end % multipareig