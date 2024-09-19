function [A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,t1,t2,t3] = ellipsoidwave_mep(n1,n2,n3,x0,y0,z0,rho,sigma,tau,opts)

%ELLIPSOIDWAVE_MEP  Discretizes ellipsoid wave system as a three-parameter eigenvalue problem
%
% [A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,t1,t2,t3] = ELLIPSOIDWAVE_MEP(n1,n2,n3,x0,y0,z0,rho,sigma,tau)
%
% Ellipsoid wave system of differential equations is
% 
% t(t-1)(t-c)F_i''(t) + 1/2*(K_2 t^2-2K_1 t+K0)F_i'(t) + (lambda -lambda0 + (mu+mu0) t + eta t^2)F_i(t)=0
%
% for i=1,2,3, where 0 < t3 < 1 < t2 < c < t1 < z0^2/b^22, where
%
% a = sqrt(z0^2-x0^2), b = sqrt(z0^2-y0^2), c = a^2/b^2, lambda0 = 1/4*((rho+tau)^2+(rho+sigma)^2*c), 
% mu0 = 1/4*(rho+sigma+tau)*(rho+sigma+tau+1), K0 = (2*rho+1)*c, K1 = (1+rho)*(1+c) + tau + sigma*c,
% K2 = 2*(rho+sigma+tau)+3
%
% Input:
%   - n1, n2, n3: number of points for X1, X2, and X3
%   - x0 < y0 < z0: semi-axis
%   - rho, sigma, tau: parameters 0 or 1 that determine one of the 8 possible types
%   - opts : options
%
% Options in opts:
%   - fp_type: numeric type to use: 'single', 'double', or 'mp' (needs MCT),
%              default is the superior type of input data (x0,y0,z0)
%
% Boundary condition on the outer boundary is Dirichlet (X1(xi0/b^2)=0)
% 
% Output:
%   - A1,B1,C1,D1 : n1 x n1 matrices for the first equation
%   - A2,B2,C2,D2 : n2 x n2 matrices for the second equation
%   - A3,B3,C3,D3 : (n3-1)x(n3-1) matrices for the third equation
%   - t1, t2, t3 : ni points for 1st, 2nd, and 3rd equation (including endpoints)
%
% See also: ELLIPSOID_EIGS, DEMO_ELLIPSOIDWAVE, DEMO_ELLIPSOIDWAVE_FIGS, BDE3MEP

% References: 
%  - M. Willatzen and L. C. Lew Yan Voon, Numerical implementation of the ellipsoidal 
%    wave equation and application to ellipsoidal quantum dots, Comput. Phys. Commun. 171 (2005) 1-18.
%  - B. Plestenjak, C. I. Gheorghiu, M. E. Hochstenbach: Spectral
%    collocation for multiparameter eigenvalue problems arising from separable
%    boundary value problems, J. Comp. Phys. 298 (2015) 585-601

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 29.11.2016: modified to be precision-independent                 
% Last revision: 29.11.2016

narginchk(9, 10);

if nargin < 10, opts = []; end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(x0,y0,z0);
end
opts.fp_type = class_t;

% Make sure all inputs are of the same numeric type.
if ~isa(x0,class_t), x0 = numeric_t(x0,class_t);  end;
if ~isa(y0,class_t), y0 = numeric_t(y0,class_t);  end;
if ~isa(z0,class_t), z0 = numeric_t(z0,class_t);  end;
if ~isa(rho,class_t), rho = numeric_t(rho,class_t);  end;
if ~isa(sigma,class_t), sigma = numeric_t(sigma,class_t);  end;
if ~isa(tau,class_t), tau = numeric_t(tau,class_t);  end;

a2 = (z0^2-x0^2);
b2 = (z0^2-y0^2);
c = a2/b2;
x10 = z0^2/b2;

lambda0 = numeric_t('1/4',class_t)*((rho+tau)^2+(rho+sigma)^2*c);
mu0 = numeric_t('1/4',class_t)*(rho+sigma+tau)*(rho+sigma+tau+1);
K0 = (2*rho+1)*c;
K1 = (1+rho)*(1+c) + tau + sigma*c;
K2 = 2*(rho+sigma+tau)+3;

p = @(x) x.*(x-1).*(x-c);
q = @(x) (K2*x.^2 - 2*K1*x + K0)/2;
r = @(x) -lambda0 + mu0*x;
s = numeric_t('-1',class_t);
t = @(x) -x;
u = @(x) -x.^2;

% Eq. 1  (c < t1 < x10)
[t1,A1,B1,C1,D1] = bde3mep(c,x10,p,q,r,s,t,u,numeric_t('[0 0;1 0]',class_t),n1,opts);
% Eq. 2  (1 < t2 < c)
[t2,A2,B2,C2,D2] = bde3mep(1,c,p,q,r,s,t,u,numeric_t('[0 0;0 0]',class_t),n2,opts);
% Eq. 3  (0 < t1 < 1)
[t3,A3,B3,C3,D3] = bde3mep(0,1,p,q,r,s,t,u,numeric_t('[0 0;0 0]',class_t),n3,opts);

