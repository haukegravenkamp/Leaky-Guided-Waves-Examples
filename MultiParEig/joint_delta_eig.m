function [lambda,report] = joint_delta_eig(Delta,opts)

%JOINT_DELTA_EIG  Solve a joined system of generalized eigenvalue problems
%
% [lambda,report] = JOINT_DELTA_EIG(Delta,opts) returns eigenvalues
% of a joint system of generalized eigenvalue problems
%
% Delta{2}*z = lambda(1)*Delta{1}*z 
% ...
% Delta{k+1}*z = lambda(k)*Delta{1}*z
%
% where Delta{1},...,Delta{k+1} are operator determinant matrices related
% to a k-parameter eigenvalue problem. If Delta{1} is nonsingular, then
% matrices inv(Delta{1})*Delta{2},...,inv(Delta{1})*Delta{k+1} commute
%
% Input:
%   - Delta : cell array of size k+1 of matrices Delta{i}, all matrices
%             have to be of the same size
%   - opts : options
%
% Options in opts:
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
%   - twosidedRQ (1): set to 0 to use one-sided Rayleigh quotients instead
%     of two-sided (only for solvers 'eig' and 'geneig') 
%   - epscluster (1e-6): relative distance between eigenvalues in a cluster,
%     set to 0 for no clustering (only for solvers 'schur' and 'genschur')
%   - rrqr (0): set to 1 to use rank revealing QR instead of SVD in the
%     staircase-type method (only for singular problems)
%   - rng_seed (0): if different from zero, rng(rng_seed) is used at start
%   - fp_type: numeric type to use ('single', 'double', or 'mp' (needs MCT),
%     use only if you want to change the default - the superior type of input data
%   - epsshape (1e-6): displays warning if shape error of simultaneous 
%     block triangularization is larger that epsshape, set to 0 for no testing
%   - all options of auxiliary functions
%
%  Outdated options that are kept for backward compatibility and now have
%  different meaning
%   - fast: no meaning - default options are used
%
% Output:
%   - lambda : matrix of size m x k, rows are eigenvalues
%   - report: details of compression method in case of a singular problem 
%       rows are [mode m n r s(1) s(r) s(r+1) choice)], 
%       where mode is 1: CR step 1, 2: CR step 2, 3: RC step 1, 4. RC step 2
%       m,n: size of Delta matrices, r: rank
%       s(1), s(r), s(r+1): singular values
%       choice: how was rank determined (see NUMRANK for details)

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% P. Holoborodko, Advanpix LLC.
% FreeBSD License, see LICENSE.txt

% BP 06.11.2022 : extracted from multipareig 
% BP 04.11.2023 : support for sparse matrices
% BP 19.12.2023 : new solvers with Rayleigh quotients
% BP 28.08.2024 : updated solvers with Rayleigh quotients

% Validate number of input parameters
narginchk(1, 2);
if nargin<2, opts=[]; end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(Delta{:});
end

% backward compatibility 
if isfield(opts,'twosideRQ') && ~isfield(opts,'twosidedRQ'), opts.twosidedRQ = opts.twosideRQ;             end

if isfield(opts,'singular'),   singular   = opts.singular;   else, singular   = 0;                         end
if isfield(opts,'epscluster'), epscluster = opts.epscluster; else, epscluster = numeric_t('1e-6',class_t); end
if isfield(opts,'rng_seed'),   rng_seed   = opts.rng_seed;   else, rng_seed   = 0;                         end
if isfield(opts,'rand_orth'),  rand_orth  = opts.rand_orth;  else, rand_orth  = 1;                         end
if isfield(opts,'solver'),     solver     = opts.solver;     else, solver = 'eig';                         end
if isfield(opts,'maxgensize'), maxgensize = opts.maxgensize; else, maxgensize = 0;                         end
if isfield(opts,'twosidedRQ'), twosidedRQ = opts.twosidedRQ; else, twosidedRQ = 1;                         end
if isfield(opts,'epsshape'),   epsshape   = opts.epsshape;   else, epsshape = 1e-6;                        end

% Default outputs
lambda = numeric_t([],class_t); 
report = numeric_t([],class_t);

k = length(Delta); % number of parameters + 1

if rng_seed
    rng(rng_seed)
end
 
for j = 1:k
    if issparse(Delta{j})
        Delta{j} = full(Delta{j});
    end
end

[m,n] = size(Delta{1});
for j = 2:k
   [m1,n1] = size(Delta{j});
   if [m1 n1] ~= [m n]
      error('Matrices must be of same size')
   end
end

% Make sure all inputs are of the same numeric type.
for j = 1:numel(Delta)
    if ~isa(Delta{j}, class_t)
         Delta{j} = numeric_t(Delta{j},class_t);
    end
end

if singular
    [Delta, report] = extract_regular_part_np(Delta, opts);
end

% Quick return in degenerate cases
if ((size(Delta{1},1) == 0) || (size(Delta{1},1) ~= size(Delta{1},2))) || ...  % no regular part was found
   ((maxgensize      > 0) && (size(Delta{1},1)  > maxgensize))                 % regular part is too large
     return
end

N = size(Delta{1},1);
out_of_shape = 0;

if rand_orth
    % we take random combination of pencils, we assume that (DeltaH,Delta{1}) 
    % has multiple eigenvalues exactly when all pencils (Delta{j},Delta{1}) 
    % have multiple eigenvalues in the same invariant subspace
    alpha = randn(k-1,1,class_t);
    alpha = alpha/norm(alpha);
    DeltaH = zeros(N,class_t);
    for j = 1:k-1
        DeltaH = DeltaH + alpha(j)*Delta{j+1};
    end
    startJ = 1;
else
    DeltaH = Delta{2};
    startJ = 2;
end

lambda = zeros(N,k-1,class_t);
if strcmp(solver,'eig')
    % we use eigenvectors of inv(Delta{1}*DeltaH) and compute eigenvalues
    % from Rayleigh quotients (standard or twosided)
    % Factorize Delta{1} only once and use factorization to solve against
    % different right-hand-sides later on. This boosts the speed slightly.
    [L,U,p] = lu(Delta{1},'vector');
    GammaH = U\(L\DeltaH(p,:));
    if twosidedRQ
        [X,D,Y] = eig(GammaH);
        leftX = Y;
    else
        [X,D] = eig(GammaH);
        leftX = X;
    end
    d0 = zeros(N,1,class_t);
    for i = 1:N
        d0(i) = leftX(:,i)'*X(:,i);
        if twosidedRQ && abs(d0(i))<1e-12
            % if left and right eigenvector are (almost) orthogoal, we switch to one-sided Rayleigh quotient for this eigenvalue
            d0(i) = 1;
            leftX(:,i) = X(:,i);
        end
    end
    if startJ == 2
        lambda(:,1) = diag(D);
    end
    for r = startJ:k-1
        % Avoid O(n^3) operations in the loop
        Gammar = (U\(L\Delta{r+1}(p,:)))*X;
        for i = 1:N
            lambda(i,r) = leftX(:,i)'*Gammar(:,i)/d0(i);
        end
    end
elseif strcmp(solver,'geneig')
    % we use eigenvectors of a GEP (DeltaH,Delta{1}) and compute eigenvalues
    % from Rayleigh quotients (standard or twosided)
    if twosidedRQ
        [X,D,Y] = eig(DeltaH, Delta{1});
        leftX = Y;
    else
        [X,D] = eig(DeltaH, Delta{1});
        leftX = X;
    end
    if startJ == 2
        lambda(:,1) = diag(D);
    end
    d0 = zeros(N,1,class_t);
    Gammar = Delta{1}*X;
    for i = 1:N
        d0(i) = leftX(:,i)'*Gammar(:,i);
        if twosidedRQ && abs(d0(i))<1e-12
            % if left and right eigenvector are (almost) orthogoal, we switch to one-sided Rayleigh quotient for this eigenvalue
            d0(i) = X(:,i)'*Gammar(:,i);
            leftX(:,i) = X(:,i);
        end
    end
    for r = startJ:k-1
        Gammar = Delta{r+1}*X;
        for i = 1:N
            lambda(i,r) = (leftX(:,i)'*Gammar(:,i))/d0(i);
        end
    end
elseif strcmp(solver,'schur') 
    % We use the Schur form of inv(Delta{1})*DeltaH (with optional reordering 
    % and clustering) to (block) triangularize other matrices. We factorize 
    % Delta{1} only once and use factorization to solve against different 
    % righ-hand-sides later on. This boosts the speed slightly. 
    [L,U,p] = lu(Delta{1},'vector');
    GammaH = U\(L\DeltaH(p,:));
    if epscluster > 0
        [Q,S,order,start,csize,tlambda] = clustered_schur(GammaH,epscluster);
    else
        [Q,S] = schur(GammaH,'complex');
        csize = 1;
        tlambda = diag(S);
    end
    if startJ == 2
        lambda(:,1) = tlambda;
    end
    if max(csize)==1
        % there are no clusters, we take diagonal elements
        for r = startJ:k-1
            % Avoid O(n^3) operations in the loop
            Gammar = (U\(L\Delta{r+1}(p,:)))*Q;
            if epsshape>0
                QG = Q'*Gammar; 
                out_of_shape(r) = norm(tril(QG,-1),'fro')/norm(QG,'fro');
                lambda(:,r) = diag(QG);
            else
                for i = 1:N
                    lambda(i,r) = Q(:,i)'*Gammar(:,i);
                end
            end
        end
    else
        % we extract eigenvalues block by block, for each block we take the
        % average eigenvalue as a representative
        for r = startJ:k-1
            mu = numeric_t([],class_t); 
            Gammar = (U\(L\Delta{r+1}(p,:)))*Q;
            S2 = Q'*Gammar;
            if epsshape>0
                out_of_shape(r) = out_of_shape_norm(S2,start,csize)/norm(S2,'fro');
            end
            for q = 1:length(start)
                if csize(q)==1
                    mu = [mu; S2(start(q),start(q))];
                else
                    partS2 = S2(start(q):(start(q)+csize(q)-1),start(q):(start(q)+csize(q)-1));
                    tmpeig = eig(partS2);
                    avgmu = sum(tmpeig)/csize(q)*ones(csize(q),1);
                    mu = [mu; avgmu];
                end
            end
            lambda(:,r) = mu;
        end
    end
elseif strcmp(solver,'genschur') 
    % We use the generalized Schur form of (DeltaH, Delta{1}) (with optional 
    % reordering and clustering) to (block) triangularize other pencils. 
    if epscluster>0
        [S1,S0,Q,Z,order,start,csize,tlambda] = clustered_qz(DeltaH,Delta{1},epscluster); %#ok<*ASGLU>
    else
        [S1,S0,Q,Z] = qz(DeltaH,Delta{1}); %#ok<*ASGLU>
        csize = 1;
        tlambda = diag(S1)./diag(S0);
    end
    if startJ == 2
        lambda(:,1) = tlambda;
    end
    if max(csize)==1
        % there are no clusters
        for r = startJ:k-1
            DZ = Delta{r+1}*Z;
            if epsshape>0
                QDZ = Q*DZ;
                out_of_shape(r) = norm(tril(QDZ,-1),'fro')/norm(QDZ,'fro');
                mu = diag(QDZ);
            else
                mu = zeros(N,1,class_t);
                 for i = 1:N
                     mu(i) = Q(i,:)*DZ(:,i);
                 end
            end
            lambda(:,r) = mu./diag(S0);
        end
    else
        % we extract eigenvalues block by block, for each block we take an
        % average eigenvalue as a representative
        for r = startJ:k-1
            mu = numeric_t([],class_t); 
            S2 = Q*Delta{r+1}*Z;
            if epsshape>0
                out_of_shape(r) = out_of_shape_norm(S2,start,csize)/norm(S2,'fro');
            end
            for p = 1:length(start)
                if csize(p)==1
                    mu = [mu; S2(start(p),start(p))/S0(start(p),start(p))];
                else
                    partS0 = S0(start(p):(start(p)+csize(p)-1),start(p):(start(p)+csize(p)-1));
                    partS2 = S2(start(p):(start(p)+csize(p)-1),start(p):(start(p)+csize(p)-1));
                    tmpeig = eig(partS2,partS0);
                    avgmu = sum(tmpeig)/csize(p)*ones(csize(p),1);
                    mu = [mu; avgmu];
                end
            end
            lambda(:,r) = mu;
        end
    end
else
    error(sprintf('Unknown solver %s, options are eig|geneig|schur|genschur',solver))
end

if epsshape>0 && max(out_of_shape)>epsshape
    warning(sprintf('Max simultaneous triangularization relative error: %0.2g',max(out_of_shape)))
end

end % joint_delta_eig

% ------------------------------------------------------------------------
% auxiliary function that returns Frobenious norm of the part outside the 
% block triangular shape
function off = out_of_shape_norm(A,start,csize)
    off = 0;
    for k=1:length(start)
        off = off + norm(A(start(k)+csize(k):end,start(k):start(k)+csize(k)-1),'fro')^2;
    end
    off = sqrt(off);
end