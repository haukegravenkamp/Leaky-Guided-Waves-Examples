function V = block_krylov_3p(B,C,F,k,A,opts)

%BLOCK_KRYLOV_3P   Orthogonal basis for generalized block Krylov subspace
% 
% V = BLOCK_KRYLOV_3P(B,C,F,k,A) returns orthogonal basis for the
% generalized block subspace, where 
% 
% k=1 : V = Lin(F, B1*F, C1*F)
% k=2 : V = Lin(F, B1*F, C1*F, B1^2*F, B1*C1*F, C1*B1*F, C1^2*F),
% k=3 : V = Lin(F, B1*F, C1*F, ..., B1^3*F, ... , C1^3*F), ..., 
%
% where B1 = inv(A)*B and C1 = inv(A)*C (if A is not supplied, B1=B, C1=C)
%
% Options in opts:
%   - svdfilter : we take only directions from SVD with relative singular values larger than svdfilter (1e-5)
%
% See also BLOCK_KRYLOV.

% References: 
%  - M. Hochstenbach, K. Meerbergen, E. Mengi, B. Plestenjak: Subspace
%    methods for 3-parameter eigenvalue problems

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Modified to be precision-independent
% Last revision: 31.1.2018

if nargin < 5, A = []; end
if nargin < 6, opts = []; end
matrixA = (nargin == 5);

if matrixA
    class_t = superiorfloat(B,C,F,A);
else
    class_t = superiorfloat(B,C,F);
end

% Make sure all inputs are of the same numeric type.
if ~isa(B,class_t), B = numeric_t(B,class_t); end;
if ~isa(C,class_t), C = numeric_t(C,class_t); end;
if ~isa(F,class_t), F = numeric_t(F,class_t); end;
if matrixA && (~isa(A,class_t)), A = numeric_t(A,class_t); end;

if isfield(opts,'svdfilter'), svdfilter = opts.svdfilter;  else svdfilter = 1e-5;  end;

mineps = numeric_t('10000*eps',class_t);

[U,tildaS,tildaV] = svd(F);
sv = diag(tildaS);
sv = sv/sv(1);
rk = sum(sv>svdfilter);
Q = U(:,1:rk);
V = Q;

zac(1) = 1;
kon(1) = rk;

for j = 1:k
    W = [B*Q C*Q];
    if matrixA
        W = A\W;
    end
    for i = 1:j
        tmp = V(:,zac(i):kon(i))'*W;
        W = W - V(:,zac(i):kon(i))*tmp;
    end
    % Repeated orhtogonalization
    for i=1:j
        tmp = V(:,zac(i):kon(i))'*W;
        W = W - V(:,zac(i):kon(i))*tmp;
    end
    % We use svd to extract only significant new directions to expand with
    [U,tildaS,tildaV] = svd(W);
    sv = diag(tildaS);
    if sv(1)<mineps
        return
    else
        sv = sv/sv(1);
        rk = sum(sv>svdfilter);
    end
    Q = U(:,1:rk);
    if rk == 0
        return
    end
    V = [V Q];
    zac(j+1) = kon(j)+1;
    kon(j+1) = kon(j)+rk;
end
    