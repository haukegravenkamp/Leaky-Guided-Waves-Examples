function [points, val, err, cand] = arma21(y,opts)

%ARMA21 returns stationary points for the ARMA(2,1) model
%
% [points, val, err, cand] = arma11(y,opts) returns real stationary points
% for the ARMA(2,1) model and objective function ||e||_2^2 
% 
% y(k) + alpha(1)*y(k-1) + alpha(2)*y(k-2) = e(k) + gamma(1)*e(k-1), k = 3,...,n
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
%   - points : matrix with real critical points [alpha(1) alpha(2) gamma(1)]
%   - val: values of the objective function 
%   - err: minimal singular values used to verify the solution
%   - cand: matrix with all eigenvalues [alpha(1) alpha(2) gamma(1)]
%
% We find critical values as eigenvalues of a rectangular two-parameter
% eigenvalue problem (A + alpha(1)*B + alpha(2)*C + gamma(1)*D + gamma(2)*E)*z = 0,
% where A,B,C,D are matrices of size (4n-5)*(4n-7)

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

% we build (4n-5)x(4n-7) matrices such that stationary points are eigenvalues 
% of the rectangular MEP (A + alfa1*B + alfa2*C + gamma1*D + gamma1^2*E) 
[A,B,C,D,E] = ARMA21_matrices(y,class_t);

% to find eigenvalues of (A + lambda*B1 + eta*B2 + mu*C + mu^2*D)x=0 
% we write the as a linear rectangular four-parameter eigenvalue problem and
% use compression

n = size(A,2);
TR = right_compression_matrix(n,3);
TL = left_compression_matrix(n+2,3);

A4 = zeros(2,class_t); A4(2,1) = 1; % A4 = [0 0;1 0];
B4 = zeros(2,class_t);              % B4 = [0 0;0 0];
C4 = zeros(2,class_t);              % C4 = [0 0;0 0];
D4 = eye(2,class_t);                % D4 = [1 0;0 1];
E4 = zeros(2,class_t); E4(1,2) = 1; % E4 = [0 1;0 0];

BCE = TL*kron3x3(B,C,E)*TR;
BCD = TL*kron3x3(B,C,D)*TR;
CDE = TL*kron3x3(C,D,E)*TR;
ACE = TL*kron3x3(A,C,E)*TR;
ACD = TL*kron3x3(A,C,D)*TR;
BDE = TL*kron3x3(B,D,E)*TR;
BAE = TL*kron3x3(B,A,E)*TR;
BAD = TL*kron3x3(B,A,D)*TR;
BCA = TL*kron3x3(B,C,A)*TR;

Delta = cell(1,5);
Delta{1} =    kron(D4,BCE) - kron(E4,BCD);                      %     BCDE
Delta{2} =   -kron(A4,CDE) - kron(D4,ACE) + kron(E4,ACD);       % -1*(ACDE)
Delta{3} =    kron(A4,BDE) - kron(D4,BAE) + kron(E4,BAD);       %     BADE
Delta{4} =   -kron(A4,BCE) + kron(E4,BCA);                      % -1* BCAE
Delta{5} =   -kron(D4,BCA) + kron(A4,BCD);                      %     BCDA

opts.singular = 1;
lambda = joint_delta_eig(Delta,opts);
alpha1 = lambda(:,1);
alpha2 = lambda(:,2);
gamma = lambda(:,3);
eta   = lambda(:,4);

msvd = zeros(length(alpha1),1); 
for k = 1:size(alpha1)
    msvd(k,:) = min(svd(A + alpha1(k)*B + alpha2(k)*C + gamma(k)*D + eta(k)*E));
end

ind = find(abs(imag(alpha1))+abs(imag(alpha2))+abs(imag(gamma))<delta);
alpha1RR = alpha1(ind);
alpha2RR = alpha2(ind);
gammaRR = gamma(ind);
sigmaRR = [];
errRR = msvd(ind);
for k = 1:size(gammaRR)
    sigmaRR(k,1) = arma21_err(y,real(alpha1RR(k)),real(alpha2RR(k)),real(gammaRR(k)),class_t);
end

points = [alpha1RR alpha2RR gammaRR];
val = sigmaRR;
cand = [alpha1 alpha2 gamma];
err = errRR;

end

function [M000,M100,M010,M001,M002] = ARMA21_matrices(y, class_t)

% Returns matrices M000, M100, M010, M001, M002 for the rectangular MEP
%   (M000 + alfa1*M100 + alfa2*M010 + gamma1*M001 + gamma1^2*M002) 
% whose eigenvalues are stationary points of the objective function for the 
% ARMA(2,1) model
 
N = length(y);
R = diag(ones(N-3,1,class_t),1) + diag(ones(N-3,1,class_t),-1);
ZB = zeros(N-2,class_t);
Id = eye(N-2,class_t);
zvec = zeros(N-2,1,class_t);
zrow = zeros(1,N-2,class_t);
y1 = y(1:N-2);
y2 = y(2:N-1);
y3 = y(3:N);

M000 = [
    y3    Id    ZB    ZB    ZB;
    y2    ZB    Id    ZB    ZB;
    y1    ZB    ZB    Id    ZB;
    zvec  R     ZB    ZB    Id; 
    0     y2'   y3'   zrow  zrow;
    0     y1'   zrow  y3'   zrow;
    0     zrow  zrow  zrow  y3'
    ];

M100 = [
    y2    ZB    ZB    ZB   ZB;
    zvec  ZB    ZB    ZB   ZB;
    zvec  ZB    ZB    ZB   ZB; 
    zvec  ZB    ZB    ZB   ZB; 
    0     zrow  y2'   zrow zrow;
    0     zrow  zrow  y2'  zrow;
    0     zrow  zrow  zrow y2'
    ];

M010 = [
    y1    ZB    ZB    ZB   ZB;
    zvec  ZB    ZB    ZB   ZB;
    zvec  ZB    ZB    ZB   ZB; 
    zvec  ZB    ZB    ZB   ZB; 
    0     zrow  y1'   zrow zrow;
    0     zrow  zrow  y1'  zrow;
    0     zrow  zrow  zrow y1'
    ];

M001 = [ 
    zvec  R     ZB    ZB   ZB;
    zvec  ZB    R     ZB   ZB;
    zvec  ZB    ZB    R    ZB;
    zvec  2*Id  ZB    ZB    R;
    0     zrow  zrow  zrow zrow;
    0     zrow  zrow  zrow zrow;
    0     zrow  zrow  zrow zrow
    ];

M002 = [ 
    zvec  Id    ZB    ZB  ZB ;
    zvec  ZB    Id    ZB  ZB ;
    zvec  ZB    ZB    Id  ZB ;
    zvec  ZB    ZB    ZB  Id ;
    0     zrow  zrow  zrow zrow;
    0     zrow  zrow  zrow zrow;
    0     zrow  zrow  zrow zrow
    ];

end

function Delta = kron3x3(A,B,C)

    Delta = kron(A,kron(B,C)-kron(C,B)) - kron(B,kron(A,C)-kron(C,A)) + kron(C,kron(A,B)-kron(B,A)); 

end

function e = arma21_err(y,alfa1,alfa2,gamma,class_t)

    N = length(y);
    TC = diag(gamma*ones(N-2,1,class_t))+diag(ones(N-3,1,class_t),1); TC(N-2,N-1) = 1;
    TA = [diag(alfa2*ones(N-2,1,class_t)) zeros(N-2,2,class_t)]+[zeros(N-2,1,class_t) diag(alfa1*ones(N-2,1,class_t)) zeros(N-2,1,class_t)] + [zeros(N-2,2,class_t) diag(ones(N-2,1,class_t))];
    err = pinv(TC)*TA*y;
    e = norm(err)^2;

end
