function [lambda,nrank,Z,d,X,Y,U,V,DA,DB,SA,SB] = singgep(A,B,opts)

%SINGGEP   Finite eigenvalues of a singular matrix pencil
%
% [lambda,nrank,Z,d,X,Y,U,V,DA,DB,SA,SB] = singgep(A,B,opts)
% computes finite regular eigenvalues of a pencil (A,B) with a projection 
% to a pencil size of a normal rank (default), rank-completing perturbation, 
% or using an augmented pencil
%
% For rank-completing perturbation or augmented pencil, if matrices A, B 
% are rectangular, we first add zero columns or rows to make them square 
% n x n matrices. 
% 
% Input
%  - A,B : matrices of the same size m x n
%
% Options in opts:
%   - nrank (-1): normal rank of A-lambda*B (if negative it is computed)
%   - show (0): display values used for the identificaton (1) or no (0)
%   - method ('project'): 'project', 'rank-complete', 'augment'
%   - tau (1e-2): perturbation size (only for rank-complete method)
%   - sep (sqrt(eps)): threshold for the separation of regular eigenvalues
%   - sepinf (100*eps): threshold for the separation of finite eigenvalues
%   - infgap (0.95): threshold for the detection of infinite eigenvalues
%   - fingap (0.01): threshold for the detection of infinite eigenvalues
%   - DA, DB: matrices of size k x k for k = max(n1,n2) - nrank for
%     prescribed eigenvalues (if not given random matrices are used) (only 
%     for rank-complete method and augment methods)
%   - SA, SB: matrices of size k x k for k = max(n1,n2) - nrank for
%     prescribed eigenvalues (only for augment method)
%   - fast_inv (0): set to 1 to compute eigenvalues of the pertrurbed 
%     pencil faster (and usually less accurate) using eig(inv(A1+B1)*A1), 
%     default is to use eig(A1,B1)
%   - fp_type: numeric type to use ('single', 'double', or 'mp' (needs MCT),
%     use only if you want to change the default - superior type of input data
%   - cplx_rnd (0): set to 1 to use complex perturbations or projections
%   - rank_est (3): number of calls rank(A+randn*B) to estimate normal rank
% 
% Output
%  - lambda: finite regular eigenvalues
%  - nrank: normal rank used in the method
%  - Z: matrix with criteria for the identification (set show=1 to see it)
%         - col. 1: index of the eigenvalue
%         - col. 2: eigenvalue of the perturbed (proj. or augm.) pencil
%         - col. 3: |s(eig)| = y'*B*x/sqrt(1+|lambda|^2)
%         - col. 4: ||V'*x||, where x is the right eigenvector
%         - col. 5: ||U'*y||, where y is the left eigenvector
%         - col. 6: gap of the eigenvalue
%         - col. 7: max(||V'*x||, ||U'*y||)
%         - col. 8: eigenvalue type: finite regular (1), infinite regular (2), 
%                   random from right singular blocks (3), random from left
%                   singular blocks (4), prescribed (5)
%  - d: eigenvalues of the perturbed pencil
%  - X, Y: right and left eigenvectors of the perturbed pencil
%  - U, V: matrices with orthonormal columns used in the method
%  - DA, DB: matrices that define prescribed eigenvalues eig(DA,DB) (only
%            for rank-complete method)
%  - SA, SB: matrices that define prescribed eigenvalues eig(SA,SB) (only
%            for augmented method)

% References: 
% 1) Algorithm 1 in M.E. Hochstenbach, C. Mehl, B. Plestenjak: Solving 
%    singular generalized eigenvalue problems by a rank-completing 
%    perturbation, SIAM J. Matrix Anal. Appl. 40 (2019) 1022-1046
% 2) Algorithms 1,2,3 in M.E. Hochstenbach, C. Mehl, B. Plestenjak: Solving 
%    singular generalized eigenvalue problems part II: projection and  
%    augmentation, arXiv 2208.01359 

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Bor Plestenjak
% 23.05.2018
% Update BP: 03.08.2022 switch to opts for setting parameters, new methods

% Tested with MPC version 4.8.8.14771

% Matlab at least 2014a is required (we need left and right eigenvectors)
if verLessThan('matlab', '8.3')
    error('Matlab at least 2014a is required')
end

narginchk(2, 3);

% Analyse user supplied options, if any
if nargin < 3, opts = []; end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(A,B);
end

if isfield(opts,'nrank'),    nrank    = opts.nrank;    else, nrank    = -1;                      end
if isfield(opts,'show'),     show     = opts.show;     else, show     = 0;                       end
if isfield(opts,'method'),   method   = opts.method;   else, method   = 'project';               end
if isfield(opts,'tau'),      tau      = opts.tau;      else, tau      = numeric_t(1e-2,class_t); end
if isfield(opts,'sep'),      sep      = opts.sep;      else, sep      = sqrt(eps(class_t));      end
if isfield(opts,'sepinf'),   sepinf   = opts.sepinf;   else, sepinf   = 100*eps(class_t);        end
if isfield(opts,'infgap'),   infgap   = opts.infgap;   else, infgap   = numeric_t(0.95,class_t); end
if isfield(opts,'fingap'),   fingap   = opts.fingap;   else, fingap   = numeric_t(1e-2,class_t); end
if isfield(opts,'DA'),       DA       = opts.DA;       else, DA       = [];                      end
if isfield(opts,'DB'),       DB       = opts.DB;       else, DB       = [];                      end
if isfield(opts,'SA'),       SA       = opts.SA;       else, SA       = [];                      end
if isfield(opts,'SB'),       SB       = opts.SB;       else, SB       = [];                      end
if isfield(opts,'fast_inv'), fast_inv = opts.fast_inv; else, fast_inv = 0;                       end
if isfield(opts,'cplx_rnd'), cplx_rnd = opts.cplx_rnd; else, cplx_rnd = 0;                       end
if isfield(opts,'scaling'),  scaling  = opts.scaling;  else, scaling  = 1;                       end
if isfield(opts,'rank_est'), rank_est = opts.rank_est; else, rank_est = 3;                       end

% Make sure all inputs are of the same numeric type.
if ~isa(A,class_t), A = numeric_t(A,class_t); end
if ~isa(B,class_t), B = numeric_t(B,class_t); end

[m,n] = size(A);

% if matrices are rectangular, we add zero columns or rows
if (strcmp(method,'rank-complete') || strcmp(method,'augment')) && m~=n
    n = max(m,n);
    A(n,n) = 0;
    B(n,n) = 0;
end

% scaling of matrices so that ||A||_1=||B||_1=1
if scaling
    alfa = norm(A,1);
    beta = norm(B,1);
else
    alfa = 1;  % no scaling so that we can compare condition numbers
    beta = 1;  % no scaling so that we can compare condition numbers
end
if (alfa == 0) || (beta == 0) % prevents NaN is one of matrices is zero
    alfa = 1;
    beta = 1;
end
A = A/alfa;
B = B/beta;
 
% if normal rank is not given, we compute it from a random linear combination of A and B 
if nrank<0 
    for r = 1:rank_est
        nr(r,1) = rank(A+randn(class_t)*B);
    end
    nrank = max(nr);
end

mn = max(m,n);
k = mn - nrank; 
if strcmp(method,'rank-complete') || strcmp(method,'augment')
    ncols = k;
else
    ncols = mn;
end
if cplx_rnd
    [U,~] = qr(randn(mn,ncols,class_t)+1i*randn(mn,ncols,class_t),0); % random complex ncols orthonormal columns
    [V,~] = qr(randn(mn,ncols,class_t)+1i*randn(mn,ncols,class_t),0); % random complex ncols orthonormal columns
else
    [U,~] = qr(randn(mn,ncols,class_t),0); % random ncols orthonormal columns
    [V,~] = qr(randn(mn,ncols,class_t),0); % random ncols orthonormal columns
    if ncols==0 % MPC incompatibility
        U = zeros(mn,0);
        V = zeros(mn,0);
    end
end
if strcmp(method, 'project') && m~=n
    if m>n
        V = V(1:n,:);
    else
        U = U(1:m,:);
    end
end
if strcmp(method,'rank-complete') || strcmp(method,'augment')
    if isempty(DA)
        DA = diag(numeric_t('1',class_t)+rand(k,1,class_t)); 
    else
        if ~isa(DA,class_t), DA = numeric_t(DA,class_t); end
    end
    if isempty(DB)
        DB = diag(numeric_t('1',class_t)+rand(k,1,class_t)); 
    else
        if ~isa(DB,class_t), DB = numeric_t(DB,class_t); end
    end
    if strcmp(method,'rank-complete') 
        % we apply rank k completing perturbation to A and B
        A1 = A + tau*U*DA*V';
        B1 = B + tau*U*DB*V';
    else % augment method
        if isempty(SA)
            SA = diag(numeric_t('1',class_t)+rand(k,1,class_t)); 
        else
            if ~isa(SA,class_t), SA = numeric_t(SA,class_t); end
        end
        if isempty(SB)
            SB = diag(numeric_t('1',class_t)+rand(k,1,class_t)); 
        else
            if ~isa(SB,class_t), SB = numeric_t(SB,class_t); end
        end
        % we augment matrices A anb B
        A1 = [A U*DA; SA*V' zeros(k,class_t)];
        B1 = [B U*DB; SB*V' zeros(k,class_t)];
    end
elseif strcmp(method,'project')
    % we project matrices A and B
    AA = U'*A*V;
    BB = U'*B*V;
    A1 = AA(k+1:end,k+1:end);
    B1 = BB(k+1:end,k+1:end);
else
    error('Unknown method, options are "rank-complete", "project", "augment"')
end
mm = size(A1,1);

% compute eigenvalues and left and right eigenvector (supported since Matlab 2014a)
if fast_inv
    % faster (not so accurate) computation using eigenvalues of inv(A1+B1)*A1
    mu = 1 + rand; % random shift (should not be zero) to make mu*A1+B1 nonsingular
    [X,D,Y] = eig((mu*A1+B1)\A1);
    d = diag(D);
    d = d./(1-mu*d);
    Y = (mu*A1'+B1')\Y;
    for i = 1:mm
        Y(:,i) = Y(:,i)/norm(Y(:,i));
    end
else    
    % slower and more reliable computation using eigenvalues of eig(A1,B1)
    [X,D,Y] = eig(A1,B1);
    d = diag(D);
    % eig(A1,B1) does not return normalized eigenvectors, thus we normalize them
    for i = 1:mm
        X(:,i) = X(:,i)/norm(X(:,i));
        Y(:,i) = Y(:,i)/norm(Y(:,i));
    end
end
B1X = B1*X;
A1X = A1*X;
s = zeros(mm,1,class_t);
a = zeros(mm,1,class_t);
b = zeros(mm,1,class_t);
e = zeros(mm,1,class_t);

for i = 1:mm
    s(i) = abs(Y(:,i)'*B1X(:,i));   % 1/condition number+
end
if strcmp(method,'rank-complete')
    VX = V'*X;
    UY = U'*Y;
elseif strcmp(method,'project')
    D = diag(d);
    U1BV2 = BB(1:k,k+1:end);
    U1AV2 = AA(1:k,k+1:end);
    U2BV1 = BB(k+1:end,1:k);
    U2AV1 = AA(k+1:end,1:k);
    MV2X = U1AV2*X-U1BV2*X*D;
    U2MX = Y'*U2AV1-D*Y'*U2BV1;
else
    VX = X(mn+1:mn+k,:);
    UY = Y(mn+1:mn+k,:);
end
if strcmp(method,'project')
    for i = 1:mm
        fac = 1+abs(d(i));
        if isinf(d(i))
            a(i) = norm(U1BV2*X(:,i));     % product V'*x
            b(i) = norm(Y(:,i)'*U2BV1);    % product U'*y
        else
            a(i) = norm(MV2X(:,i))/fac;    % product V'*x
            b(i) = norm(U2MX(i,:))/fac;    % product U'*y
        end
        e(i,1) = max(a(i,1),b(i,1));       % criterion for separation
    end
else    
    for i = 1:mm
        a(i) = norm(VX(:,i));       % product V'*x
        b(i) = norm(UY(:,i));       % product U'*y
        e(i) = max(a(i),b(i));      % criterion for separation
    end
end

s = s./sqrt(1+abs(d).^2);
d = d*alfa/beta; % eigenvalues (scaled back to original (A,B))
eigv_type = zeros(mm,1); 

cand_regular = e<sep;             % criterion for regular eigenvalues 
pos_regular = find(cand_regular); % positions of regular eigenvalues
pos_fake = find(1-cand_regular);  % positions of fake eigenvalues

eigv_type(pos_regular) = 1;
for j = pos_fake(:)'
    if a(j)/b(j)<sep            
        % if ||V'*x||=0 and ||U'*y|| nonzero this is a random eigenvalue from the right singular block
        eigv_type(j) = 3;
    elseif b(j)/a(j)<sep
        % if ||U'*y||=0 and ||V'*x|| nonzero this is a random eigenvalue from the left singular block
        eigv_type(j) = 4;
    else
        % if both ||V'*x|| and ||U'*y|| are nonzero, this is a prescribed eigenvalue
        eigv_type(j) = 5;
    end
end

gap = reldist(d);       % eigenvalue gaps

cand_infinite = (s<sep).*(gap>infgap) + (s<=sepinf).*(gap>fingap); % criteria for infinite eigenvalues
ind_inf = find(cand_regular.*cand_infinite); % positions of infinite regular eigenvalues 
eigv_type(ind_inf) = 2;

[~, ord1] = sort(-s); % we sort first by s (so that regular eigenvalues will be sorted by s)
[eigv_type, ord2] = sort(eigv_type(ord1)); % next we sort by eigenvalue type
finord = ord1(ord2);
a = a(finord);
b = b(finord);
d = d(finord);
s = s(finord);
gap = gap(finord);
e = e(finord);
X = X(:,finord);
Y = Y(:,finord);
lambda = d(eigv_type==1);          % finite regular eigenvalues

Z = [(1:mm)' d s a b gap e eigv_type]; 

if show
    disp(' ')
    disp(' i  |     re(lambda)     im(lambda) |     s(eig) |   ||V''*x|| |   ||U''*y|| |    gap     | max(Vx,Uy) | type')
    disp('--------------------------------------------------------------------------------------------------------------')
    for i=1:mm
        j = i;
        fprintf('%3d | %14.6e %14.6e | %10.2e | %10.2e | %10.2e | %10.2e | %10.2e |  %d\n',i,real(d(j)),imag(d(j)),s(j),a(j),b(j),gap(j), e(j),eigv_type(j));
    end
end

end % singgep

% auxiliary function that computes gaps of the eigenvalue
function gap = reldist(x)
    n = length(x);
    gap = [];
    if n>0
        M = x*ones(1,n)-ones(n,1)*(x.')+diag(Inf*ones(n,1));
        d = min(abs(M));
        a = d.';         % absolute gaps
        y = a./abs(x);   % relative gaps 
        gap = min(a,y);
        gap = a./sqrt(1+abs(x).^2);
        gap(isnan(gap)) = 1;
    end
end