function [A1,B1,C1,A2,B2,C2,z1,z2] = hydrogen_mep(n1,n2,b,R,k,opts)

%HYDROGEN_MEP  Discretizes system of Schroedinger equation for hydrogen
%molecular ion in 2D as a two-parameter eigenvalue problem
%
% [A1,B1,C1,A2,B2,C2,z1,z2] = HYDROGEN_MEP(n1,n2,b,R,k)
% transforms system of differential equations 
% 
% (u^2-1) f''(u) + (2k+1)u f'(u) + (k + 2Ru + 1/2*R^2*mu*(u^2-1) - lambda) f(u) = 0,  1 < u < Inf,
% (v^2-1) g''(v) + (2k+1)v g'(v) + (k + 1/2*R^2*mu*(v^2-1) - lambda) g(v) = 0,  -1 < v < 1,
%
% into a two-parameter eigenvalue problem using Laguerre and Chebyshev collocation
%
% Input:
%   - n1, n2: number of points for the first and the second equation
%   - b: scaling parameter for Laguerre collocation
%   - R, k : parameter for the equation (k is 0 or 1)
%   - opts : options
%
% Options in opts:
%   - fp_type: numeric type to use: 'single', 'double', or 'mp' (needs MCT),
%              default is the superior type of input data (x0,y0,z0)
%
% Output:
%   - A1,B1,C1 : n1 x n1 matrices for the first equation
%   - A2,B2,C2 : n2 x n2 matrices for the second equation
%   - z1, z2 : ni points for 1st and 2nd equation (including endpoints)
%
% See also: HYDROGEN_MODES, DEMO_HYDROGEN, BDE2MEP

% Package dmsuite by J.A.C Weideman is required (or other chebdif and lagdif methods)

% References: 
%  - Patil, Hydrogen molecular ion and molecule in two dimensions, 
%    Journal of Chemical Physics 118, (2003) 2197-2205.
%  - B. Plestenjak, C. I. Gheorghiu, M. E. Hochstenbach: Spectral
%    collocation for multiparameter eigenvalue problems arising from separable
%    boundary value problems, J. Comp. Phys. 298 (2015) 585-601

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 03.12.2016: modified to be precision-independent                 
% Last revision: 03.12.2016

narginchk(5, 6);

if nargin < 6, opts = []; end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(b,R);
end
opts.fp_type = class_t;

% Make sure all inputs are of the same numeric type.
if ~isa(b,class_t), b = numeric_t(b,class_t);  end;
if ~isa(R,class_t), R = numeric_t(R,class_t);  end;

% radial equation
% -----------------------------------------------------------
% radial equation F(u) = f(u)(u^2-1)^(m/2) on [1,Inf]
% boundary conditions: (2m+1)*f'(1) = -(k + 2R - Lambda)f(1) and f(Inf) = 0

[x1,tmpD] = lagdif_mp(n1,2,b,class_t); % Laguerre points
z1 = x1 + 1; % shift to [1, Infty]
D0 = eye(n1,class_t);
D1 = tmpD(:,:,1);
D2 = tmpD(:,:,2);

diagP1 = diag(z1.^2-1); 
diagQ1 = diag((2*k+1)*z1); 
diagR1 = diag(2*R*z1 + k); 
diagS1 = eye(n1,class_t); 
diagT1 = diag(-1/2*R^2*(z1.^2-1)); 

A1 = diagP1*D2 + diagQ1*D1 + diagR1;
B1 = diagS1;
C1 = diagT1;

% boundary condition in point 1 (this is first row): 
A1(1,:) = (2*k+1)*D1(1,:) + (2*R + k)*D0(1,:);
B1(1,:) = D0(1,:);
C1(1,:) = zeros(1,n1,class_t);

% angular equation
% -----------------------------------------------------------
% angular equation G(y) = (1-y^2)^(m/2)g(y) on [-1,1]
% boundary conditions:  (2m+1)g'(1) + (m - Lambda) g(1) = 0
%                      -(2m+1)g'(-1) +(m - lambda) g(-1) = 0

a2 = numeric_t('-1',class_t); 
b2 = numeric_t('1',class_t);
p2 = @(y) y.^2-1;
q2 = @(y) (2*k+1)*y;
r2 = k;
s2 = numeric_t('1',class_t);
t2 = @(y) -1/2*R^2*(y.^2-1);
[z2,A2,B2,C2,G2,k2,r2] = bde2mep(a2,b2,p2,q2,r2,s2,t2,numeric_t('[0 0; 0 0]',class_t),n2);
