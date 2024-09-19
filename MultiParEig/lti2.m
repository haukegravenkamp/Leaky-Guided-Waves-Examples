function [points, val, err, cand] = lti2(y,opts)

%LLTI2 returns stationary points for the LTI(2) model
%
% [points, val, err, cand] = lti2(y,opts) returns real stationary points
% for the LTI(2) model and objective function ||e||_2^2, where we are
% looking for parameters alpha(1) and alpha(2) for teh best 2-norm
% approximation of y by z, such that its elements satisfy the difference
% equation 
%
% z(k+2) + alpha(1)*z(k+1) + alpha(2)*z(k) = 0, k = 1,...,n-2
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
%   - points : matrix with real critical points [alpha(1) alpha(2)]
%   - val: values of the objective function 
%   - err: minimal singular values used to verify the solution
%   - cand: matrix with all eigenvalues [alpha(1) gamma(1)]
%
% We find critical values as eigenvalues of a rectangular two-parameter
% eigenvalue problem (A00 + alpha(1)*A10 + alpha(2)*A01 + alpha(1)^2*A20 + alpha(1)*alpha(2)*A11 + alpha(2)^2*A20)*z = 0,
% where A00,A10,A01,A20,A11 (A02=A20) are matrices of size (3n-4)*(3n-5)

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

% we build (3n-4)x(3n-5) matrices such that stationary points are eigenvalues 
% rectangular MEP (M00 + alpha(1)*M10 + alpha(2)*M01 + alpha(1)^2*M20 + alpha(1)*alpha(2)*M11 + alpha(2)^2*M20)*z=0 
[M00,M10,M01,M20,M11,M02] = LTI2_matrices(y);

% linearization for alpha_3 = alpha_1*alpha2
A3 = zeros(2,class_t); A3(2,1) = 1; % A3 = [0 0;1 0];
B3 = zeros(2,class_t); B3(1,1) = 1; % B3 = [1 0; 0 0];
C3 = zeros(2,class_t); C3(2,2) = 1; % C3 = [0 0; 0 1];
D3 = zeros(2,class_t); D3(1,2) = 1; % D3 = [0 1; 0 0];
E3 = zeros(2,class_t);

% linearization for alpha_4 = alpha_1^2 + alpha2^2
A4 = zeros(2,class_t); A4(2,1) = 1; % A4 = [0 0; 1 0];
B4 = eye(2,class_t);                % B4 = [1 0; 0 1];
C4 = eye(2,class_t);                % C4 = [1 0; 0 1];
D4 = zeros(2,class_t); D4(1,2) = 2; % D4 = [0 2; 0 0];
E4 = zeros(2,class_t); E4(1,2) = 1; % E4 = [0 1; 0 0];

W = cell(4,5);
W{1,1} = -M00; W{1,2} = M10; W{1,3} = M01; W{1,4} = M11; W{1,5} = M02; 
W{2,1} = -M00; W{2,2} = M10; W{2,3} = M01; W{2,4} = M11; W{2,5} = M02; 
W{3,1} = -A3; W{3,2} = B3; W{3,3} = C3; W{3,4} = D3; W{3,5} = E3;
W{4,1} = -A4; W{4,2} = B4; W{4,3} = C4; W{4,4} = D4; W{4,5} = E4;

RectDelta = multipar_delta(W);

n = size(M00,2);
TR = right_compression_matrix(n,2);
TL = left_compression_matrix(n+1,2);
bigTL = kron(TL,eye(4));
bigTR = kron(TR,eye(4));

Delta = cell(1,5);
for k = 1:5
    Delta{k} = bigTL*RectDelta{k}*bigTR;
end

opts.singular = 1;
lambda = joint_delta_eig(Delta,opts);
alpha1 = lambda(:,1);
alpha2 = lambda(:,2);

msvd = [];
for k = 1:length(alpha1)
    msvd(k,1) = min(svd(M00 + alpha1(k)*M10 + alpha2(k)*M01 + alpha1(k)^2*M20 + alpha1(k)*alpha2(k)*M11 + alpha2(k)^2*M02));
end

ind = find(abs(imag(alpha1))+abs(imag(alpha2))<1e-6);
alpha1RR = alpha1(ind);
alpha2RR = alpha2(ind);
sigmaRR = [];
errRR = msvd(ind);
for k = 1:size(alpha2RR)
    sigmaRR(k,1) = LTI2_err(y,real(alpha1RR(k)),real(alpha2RR(k)),class_t);
end

res = [alpha1RR alpha2RR errRR sigmaRR];
points = [alpha1RR alpha2RR];
val = sigmaRR;
err = errRR;
cand = [alpha1 alpha2];

end

function [M00,M10,M01,M20,M11,M02] = LTI2_matrices(y)

% [M00,M10,M01,M20,M11,M02] = LTI_matrices(y) returns (3N-4)x(3N-5) 
% matrices M00, M10, M01, M20, M11, M02 for the rectangular MEP
% (M00 + l*M10 + u*M01 + l^2*M20 + l*u*M11 + u^2*M02)z = 0
% for the optimal LTI(2) identification problem

class_t = superiorfloat(y);

N = length(y);
y = y(:);
y1 = y(1:N-2);
y2 = y(2:N-1);
y3 = y(3:N);
Id = eye(N-2,class_t);
ZM = zeros(N-2,class_t);
Zr = zeros(1,N-2,class_t);
Zc = zeros(N-2,1,class_t);

P010 = diag(ones(N-3,1,class_t),1) + diag(ones(N-3,1,class_t),-1);
if N==3
    P100 = numeric_t(0,class_t);
else
    P100 = diag(ones(N-4,1,class_t),2) + diag(ones(N-4,1,class_t),-2);
end

M00 = [y3   Id     ZM    ZM
       y2   P010   Id    ZM
       y1   P100   ZM    Id
       0    y2.'   y3.'  Zr
       0    y1.'   Zr    y3.'];

M10 = [y2   P010   ZM    ZM
       Zc   2*Id   P010  ZM
       Zc   P010   ZM    P010
       0    Zr     y2.'  Zr
       0    Zr     Zr    y2.'];

M01 = [y1   P100   ZM    ZM
       Zc   P010   P100  ZM
       Zc   2*Id   ZM    P100
       0    Zr     y1.'  Zr
       0    Zr     Zr    y1.'];

M11 = [Zc   P010   ZM    ZM
       Zc   ZM     P010  ZM
       Zc   ZM     ZM    P010
       0    Zr     Zr    Zr
       0    Zr     Zr    Zr];

M20 = [Zc   Id     ZM    ZM
       Zc   ZM     Id    ZM
       Zc   ZM     ZM    Id
       0    Zr     Zr    Zr
       0    Zr     Zr    Zr];

M02 = M20;

end

function e = LTI2_err(y,alpha1,alpha2,class_t)

    N = length(y);

    tmp1 = [alpha2*eye(N-2,class_t) zeros(N-2,2,class_t)];
    tmp2 = [zeros(N-2,1,class_t) alpha1*eye(N-2,class_t) zeros(N-2,1,class_t)];
    tmp3 = [zeros(N-2,2,class_t) eye(N-2,class_t)];
    TA = tmp1 + tmp2 + tmp3;
    DC = TA*TA';

    err = TA'*(DC\(TA*y));
    e = norm(err)^2;

end