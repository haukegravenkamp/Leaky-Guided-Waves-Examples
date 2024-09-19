function [lambda,mu,X1,X2,flag] = twopareigs_ira(A1,B1,C1,A2,B2,C2,k,opts)

%TWOPAREIGS_IRA   Implicitly restarted inverse Arnoldi for two-parameter eigenvalue problem
%
% [lambda,mu,X1,X2] = TWOPAREIGS_IRA(A1,B1,C1,A2,B2,C2,k,opts)
% returns k eigenvalues with the smallest |mu| of the two-parameter
% eigenvalue problem
%
% A1*x = lambda*B1*x + mu*C1*x
% A2*y = lambda*B2*y + mu*C2*y
%
% using implicit restarted inverse Arnoldi in Matlab routine eigs
% on the generalized eigenvalue problem 
% Delta2*z = lambda*Delta0*z, where z = kron(x,y) and
% Delta0 = kron(B1,C2) - kron(C1,B2)
% Delta2 = kron(B1,A2) - kron(A1,B2)
% 
% When building Krylov subspace, Bartels-Stewart method is used to solve 
% the Sylvester equation related to linear system Delta2*w = Delta0*z
% 
% Input:
%   - A1, B1, C1, A2, B2, C2 : real matrices
%   - k : number of eigenvalues(6)
%   - opts: options (see below)
%
% Output: 
%   - lambda, mu: eigenvalue parts (eigenvalues are (lambda(j),mu(j))
%   - X1, X2: components of decomposable right eigenvectors (eigenvector is kron(X1(:,j),X2(:,j)) such that 
%       (A1-lambda(j)*B1-mu(j)*C1)X1(:,j)=0
%       (A2-lambda(j)*B2-mu(j)*C2)X2(:,j)=0
%   - flag: convergence (0), no convergence (1)
%
% Options in opts:
%   - divA : (default 1) : divide by A1 and A2, or by B1 and B2 (0)
%   - refine : number of TRQ refinement steps in the end (1)
%   - epscluster : distance between eigenvalues in the same cluster (1e-4)%   
%   - init_double: (default 1) if the numerical type used is not double,
%     problem is solved first in double precision and then refined, if you
%     set to 0, then iterative methods will use multiprecision, which is
%     much slower
%   - refine_mp : TRQ refinement steps in multiprecision (10)
%   - all options for the eigs method (except isreal)
%
% For Matlab below 2014a package lapack from MatlabCentral is required 
% for faster evaluation 
%
% See also: TWOPAREIG, TWOPAREIGS, TWOPAREIGS_KS, TWOPAREIGS_JD,
% TWOPAREIGS_SI, THREEPAREIGS.

% Reference: K. Meerbergen, B. Plestenjak: A Sylvester-Arnoldi type method 
% for the generalized eigenvalue problem with two-by-two operator determinants, 
% Numer. Linear Algebra Appl. 22 (2015) 1131-1146

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 02.02.2018: modified for MCT
% BP 30.03.2015: bug with complex matrices solved
% BP 25.08.2015: in Matlab 2014a or later the builtin function sylvester is
% used, renamed from twopareigs, which is now the main driver 

% Last revision: 2.2.2018

narginchk(6, 8);

if nargin<7, k = 6; end
if nargin<8, opts=[]; end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(A1,B1,C1,A2,B2,C2);
end
if isfield(opts,'init_double'), init_double = opts.init_double;  else init_double = 1; end
% class_ks is the class that is used in the procedure
if init_double
    class_ks = 'double';
else
    class_ks = class_t;
end

if isfield(opts,'divA'),       divA = opts.divA;              else divA = 1;          end
if isfield(opts,'refine'),     refine = opts.refine;          else refine = 1;        end
if isfield(opts,'epscluster'), epscluster = opts.epscluster;  else epscluster = 1e-4; end
if isfield(opts,'refine_mp'),  refine_mp = opts.refine_mp;    else refine_mp = 10;    end

% turn off diagnostic information level in eigs in earlier versions of Matlab
if ~isfield(opts,'disp'),      opts.disp=0;  end; 

n1 = size(A1,1);
n2 = size(A2,1);

% for 1x1 matrices we do not need eigs
if n1*n2==1
    [lambda,mu,X1,X2] = twopareig(A1,B1,C1,A2,B2,C2);
    return
end

if isreal(A1) && isreal(B1) && isreal(C1) && isreal(A2) && isreal(B2) && isreal(C2)
    realMEP = 1; 
else
    realMEP = 0;
    opts.isreal = false;
end

if ~verLessThan('matlab', '8.3') || strcmpi(class_ks,'mp')
    uselapack = 0; % Matlab at least 2014a or MCT
elseif exist('lapack','file') 
    uselapack = 1; % old Matlab with lapack package
else
    uselapack = 2; % old Matlab without lapack package
end

if init_double
    MA1 = double(A1); MB1 = double(B1); MC1 = double(C1);
    MA2 = double(A2); MB2 = double(B2); MC2 = double(C2);
else
    MA1 = A1; MB1 = B1; MC1 = C1; 
    MA2 = A2; MB2 = B2; MC2 = C2;
end

% Schur decomposition for Bartels-Stewart
if uselapack<2
    % Schur decomposition (real or complex) for Bartels-Stewart
    if divA
        [U1,R1] = schur(MA2\MB2);
        [U2,R2] = schur(-transpose(MB1)/transpose(MA1));
    else
        [U1,R1] = schur(MB2\MA2);
        [U2,R2] = schur(-transpose(MA1)/transpose(MB1));
    end
else
    % in old Matlab we use complex Schur even for real matrices to solve the Sylvester equaiton
    if divA
        [U1,R1] = schur(MA2\MB2,'complex');
        [U2,R2] = schur(-transpose(MB1)/transpose(MA1),'complex');
    else
        [U1,R1] = schur(MB2\MA2,'complex');
        [U2,R2] = schur(-transpose(MA1)/transpose(MB1),'complex');
    end
end

if ~realMEP
    opts.isreal = 0;
end

if strcmpi(class_ks,'mp')
    [Z,D,flag] = mpeigs(@multGamma2,n1*n2,k,'SM',opts);
else
    [Z,D,flag] = eigs(@multGamma2,n1*n2,k,'SM',opts);
end

% -- Function for multiplication by inv(Delta2)*Delta0 in complex case
function y = multGamma2(x)
    W = reshape(x,n2,n1);
    FF = MC2*W*transpose(MB1) - MB2*W*transpose(MC1); % f = Delta0*w
    if divA
        y = reshape(sylv(U1,R1,U2,R2, -MA2\FF/transpose(MA1),uselapack),n1*n2,1); % wt = Delta2\Delta0*w
    else
        y = reshape(sylv(U1,R1,U2,R2, MB2\FF/transpose(MB1),uselapack),n1*n2,1); % wt = Delta2\Delta0*w
    end
    if realMEP && (uselapack==2) % we take only the real part if a real system was solved using complex Schur 
        y = real(y);
    end
end

% we have to cluster eigenvalues in case we have multiple mu
[order,clstart,clsize,mu] = clusters(diag(D),epscluster);
Z = Z(:,order); % eigenvectors sorted in clusters
for c = 1:length(clstart)
   % for each cluster we compute matrices Gi = W'*Deltai*W
   % that give as lambda and correct eigenvectors
   % in cluster of size 1 this is just Rayleigh quotient
   j1 = clstart(c);
   j2 = clstart(c)+clsize(c)-1;
   [G0,G1,tilda] = mult_delta(MA1,MB1,MC1,MA2,MB2,MC2,Z(:,j1:j2),Z(:,j1:j2));
   [tmpX,tmpD] = eig(G1,G0);
   Z(:,j1:j2) = Z(:,j1:j2)*tmpX;
   lambda(j1:j2,1) = diag(tmpD);
end

k = length(mu);
for j=1:k
    if refine
       mz = reshape(Z(:,j),n2,n1);
       [tmpU,tilda,tmpV] = svd(mz);
       refx1 = conj(tmpV(:,1));
       refx2 = tmpU(:,1);
       [refl,refu,refx1,refx2] = trqi(MA1,MB1,MC1,MA2,MB2,MC2,refx1,refx2,refine,eps);
       mu(j,1) = refu;
       lambda(j,1) = refl;
       X1(:,j) = refx1;
       X2(:,j) = refx2;
    else    
       % extraction of eigenvectors 
       X1(:,j) = min_sing_vec(MA1-lambda(j)*MB1-mu(j)*MC1);
       X2(:,j) = min_sing_vec(MA2-lambda(j)*MB2-mu(j)*MC2);
    end
end

if init_double && ~strcmpi(class_t,'double')
    % TRQI refinement in higher precision
    RX1 = zeros(n1,k,class_t);
    RX2 = zeros(n2,k,class_t);
    Rlambda = zeros(k,1,class_t);
    Rmu = zeros(k,1,class_t);
    for j = 1:length(lambda)
        [Rlambda(j,1),Rmu(j,1),RX1(:,j),RX2(:,j)] = trqi(A1,B1,C1,A2,B2,C2,numeric_t(X1(:,j),class_t),numeric_t(X2(:,j),class_t),refine_mp,numeric_t('eps',class_t));
    end
    X1 = RX1;
    X2 = RX2;
    lambda = Rlambda;
    mu = Rmu;
end

end