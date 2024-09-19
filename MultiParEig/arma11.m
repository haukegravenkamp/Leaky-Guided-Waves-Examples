function [points, val, err, cand] = arma11(y,opts)

%ARMA11 returns stationary points for the ARMA(1,1) model
%
% [points, val, err, cand] = arma11(y,opts) returns real stationary points
% for the ARMA(1,1) model and objective function ||e||_2^2 
% 
% y(k) + alpha(1)*y(k-1) = e(k) + gamma(1)*e(k-1), k = 2,...,n
%
% for a given vector y of length n.
%
% Input:
%   - y : real vector of size n
%   - opts : options 
%
% Options in opts:
%   - fp_type: numeric type to use ('single', 'double', or 'mp' (needs MCT),
%     use only if you want to change the default - the superior type of input data
%   - delta (sqrt(eps)): treshold for the extraction of real solutions
%   - all options for joint_delta_eig
%
% Output:
%   - points : matrix with real critical points [alpha(1) gamma(1)]
%   - val: values of the objective function 
%   - err: minimal singular values used to verify the solution
%   - cand: matrix with all eigenvalues [alpha(1) gamma(1)]
%
% We find critical values as eigenvalues of a rectangular two-parameter
% eigenvalue problem (A + alpha(1)*B + gamma(1)*C + gamma(2)*D)*z = 0,
% where A,B,C,D are matrices of size (3n-1)*(3n-2)

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% P. Holoborodko, Advanpix LLC.
% FreeBSD License, see LICENSE.txt

% Bor Plestenjak 2022

% See: M.E.Hochstenbach, T.Kosir, B.Plestenjak: Numerical methods for rectangular 
% multiparameter eigenvalue problems, with applications to finding optimal 
% ARMA and LTI models. Numer Linear Algebra Appl. 2023; e2540

narginchk(1,2);

% Analyse user supplied options, if any
if nargin < 2, opts = []; end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(y);
end

if isfield(opts,'delta'),    delta   = opts.delta;    else,   delta = sqrt(eps(class_t));   end
if ~isfield(opts,'epscluster'), opts.epscluster = 1e-3;  end

y = y(:);
if ~isa(y,class_t), y = numeric_t(y,class_t); end

% we build (3n-1)x(3n-2) matrices such that stationary points are eigenvalues 
% of the rectangular MEP (A + alfa1*B + gamma1*C + gamma1^2*D) 
[A,B,C,D] = ARMA11_matrices(y,class_t);

% to find eigenvalues of (A + lambda*B + mu*C + mu^2*D)x=0 we write the
% problem as a linear rectangular three-parameter eigenvalue problem and
% use compression

n = size(A,2);
TR = right_compression_matrix(n,2);
TL = left_compression_matrix(n+1,2);

A3 = zeros(2,class_t); A3(2,1) = 1; % A3 = [0 0;1 0];
B3 = zeros(2,class_t);              % B3 = [0 0;0 0];
C3 = eye(2,class_t);                % C3 = [1 0;0 1];
D3 = zeros(2,class_t); D3(1,2) = 1; % D3 = [0 1;0 0];

BC = TL*(kron(B,C)-kron(C,B))*TR;
BD = TL*(kron(B,D)-kron(D,B))*TR;
CD = TL*(kron(C,D)-kron(D,C))*TR;
AC = TL*(kron(A,C)-kron(C,A))*TR;
AD = TL*(kron(A,D)-kron(D,A))*TR;
BA = TL*(kron(B,A)-kron(A,B))*TR;

Delta = cell(1,4);
Delta{1} = -(kron(BC,D3) - kron(BD,C3));
Delta{2} = kron(AC,D3) - kron(AD,C3) + kron(CD,A3);
Delta{3} = kron(BA,D3) - kron(BD,A3) + kron(AD,B3);
Delta{4} = kron(BC,A3) - kron(BA,C3) + kron(CD,B3);

opts.singular = 1;
lambda = joint_delta_eig(Delta,opts);
alpha = lambda(:,1);
gamma = lambda(:,2);
eta = lambda(:,3);

msvd = zeros(length(alpha),1); 
for k = 1:length(alpha)
    msvd(k,:) = min(svd(A+alpha(k)*B+gamma(k)*C+eta(k)*D));
end

ind = find(abs(imag(alpha))+abs(imag(gamma))<delta);
alphaRR = alpha(ind);
gammaRR = gamma(ind);
msvdRR = msvd(ind);
sigmaRR = zeros(length(gammaRR),1);
for k = 1:length(gammaRR)
    sigmaRR(k,1) = arma11_err(y,real(alphaRR(k)),real(gammaRR(k)),class_t);
end

points = [alphaRR gammaRR];
val = sigmaRR;
cand = [alpha gamma];
err = msvdRR;

end

function [M00,M10,M01,M02] = ARMA11_matrices(y,class_t)

% Returns matrices M00, M10, M01, M11 for the rectangular MEP
%   (M00 + alfa1*M10 + gamma1*M01 + gamma1^2*M02) 
% whose eigenvalues are stationary points of the objective function for the 
% ARMA(1,1) model
 
N = length(y);
R = diag(ones(N-2,1,class_t),1) + diag(ones(N-2,1,class_t),-1);
ZB = zeros(N-1,class_t);
Id = eye(N-1,class_t);
zvec = zeros(N-1,1,class_t);
zrow = zeros(1,N-1,class_t);
y1 = y(1:N-1);
y2 = y(2:N);

M00 = [
    y2    Id    ZB    ZB;
    y1    ZB    Id    ZB;
    zvec  R     ZB    Id; 
    0     y1'   y2'   zrow;
    0     zrow  zrow  y2'
    ];

M10 = [
    y1    ZB    ZB    ZB;
    zvec  ZB    ZB    ZB;
    zvec  ZB    ZB    ZB; 
    0     zrow  y1'   zrow;
    0     zrow  zrow  y1'
    ];

M01 = [ 
    zvec  R     ZB    ZB;
    zvec  ZB    R     ZB;
    zvec  2*Id  ZB    R;
    0     zrow  zrow  zrow;
    0     zrow  zrow  zrow
    ];

M02 = [ 
    zvec  Id    ZB    ZB;
    zvec  ZB    Id    ZB;
    zvec  ZB    ZB    Id;
    0     zrow  zrow  zrow;
    0     zrow  zrow  zrow
    ];

end

function e = arma11_err(y,alpha,gamma,class_t)
% returns value of the objective function for the ARMA(1,1) model

    N = length(y);
    TC = diag(gamma*ones(N-1,1,class_t))+diag(ones(N-2,1,class_t),1); TC(N-1,N) = 1;
    TA = diag(alpha*ones(N-1,1,class_t))+diag(ones(N-2,1,class_t),1); TA(N-1,N) = 1;
    err = pinv(TC)*TA*y;
    e = norm(err)^2;

end
