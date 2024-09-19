function [A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,t1,t2,t3] = baer_mep(n1,n2,n3,c,b,c0,b0,opts)

%BAER_MEP  Discretizes system of Baer wave equations as a 3-parameter eigenvalue problem
%
% [A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,t1,t2,t3] = baer_mep(n1,n2,n3,c,b,nu0,mu0)
%
% Baer wave system of differential equations is
% 
% (ti-b)(ti-c)F_i''(ti) + 1/2*(2*ti-(b+c))F_i'(ti) + (lambda + mu*ti + eta*ti^2)F_i(ti)=0
%
% for i=1,2,3, where c0 < t1 < c < t2 < b < t3 < b0
%
% System corresponds to separated Helmholtz equation in (confocal) paraboloidal
% coordinates (t1,t2,t3), where t1 < c < t2 < b < t3. If [lambda mu eta] 
% is an eigenvalues, then eta = omega^2, where omega is the wavenumber from 
% the Helmholtz equation.
% 
% We get eigenmodes for the body bounded by two elliptic paraboloids 
% xi1 = c0 and xi3 = b0 and Dirichlet (zero) boundary condition
%
% Input:
%   - n1, n2, n3: number of points for t1, t2, and t3
%   - c < b: free parameters of the paraboloidal coordinate system
%   - c0 < c and b < b0: intersection with z plane of the elliptic paraboloids
%   - opts : options
%
% Options in opts:
%   - fp_type: numeric type to use: 'single', 'double', or 'mp' (needs MCT),
%              default is the superior type of input data (c,b,c0,b0)
%
% Boundary condition on the outer boundary is Dirichlet (F1(c0)=0 and F3(b0)=0)
% 
% Output:
%   - A1,B1,C1,D1 : (n1-1) x (n1-1) matrices for the first equation
%   - A2,B2,C2,D2 : n2 x n2 matrices for the second equation
%   - A3,B3,C3,D3 : (n3-1)x(n3-1) matrices for the third equation
%   - t1, t2, t3 : ni points for 1st, 2nd, and 3rd equation (including endpoints)

% References: 
%  - M. Hochstenbach, K. Meerbergen, E. Mengi, B. Plestenjak: Subspace
%    methods for 3-parameter eigenvalue problems, arXiv 1802:07386

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Supports MCT

% Bor Plestenjak 26.02.2018

narginchk(7, 8);

if nargin<8; 
    opts = [];
end

if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(c,b,c0,b0);
end
opts.fp_type = class_t;

% Make sure all inputs are of the same numeric type.
if ~isa(c,class_t),  c  = numeric_t(c,class_t);  end;
if ~isa(b,class_t),  b  = numeric_t(b,class_t);  end;
if ~isa(c0,class_t), c0 = numeric_t(c0,class_t);  end;
if ~isa(b0,class_t), b0 = numeric_t(b0,class_t);  end;

% Eq. 1  F1(t1) c0 < t1 < c
p1 = @(x) (b-x).*(c-x);
q1 = @(x) 1/2*(2*x-(b+c));
r1 = numeric_t('0',class_t);
s1 = numeric_t('-1',class_t);
t1 = @(x) -x;
u1 = @(x) -x.^2;
bc1 = numeric_t('[1 0; 0 0]',class_t); 
[t1,A1,B1,C1,D1] = bde3mep(c0,c,p1,q1,r1,s1,t1,u1,bc1,n1,opts);

% Eq. 2  F2(t2) c < t2 < b
p2 = @(x) (b-x).*(c-x);
q2 = @(x) 1/2*(2*x-(b+c));
r2 = numeric_t('0',class_t);
s2 = numeric_t('-1',class_t);
t2 = @(x) -x;
u2 = @(x) -x.^2;
bc2 = numeric_t('[0 0; 0 0]',class_t); 
[t2,A2,B2,C2,D2] = bde3mep(c,b,p2,q2,r2,s2,t2,u2,bc2,n2,opts);

% Eq. 3 F3(t3) b < t3 < b0
p3 = @(x) (b-x).*(c-x);
q3 = @(x) 1/2*(2*x-(b+c));
r3 = numeric_t('0',class_t);
s3 = numeric_t('-1',class_t);
t3 = @(x) -x;
u3 = @(x) -x.^2;
bc3 = numeric_t('[0 0; 1 0]',class_t); 
[t3,A3,B3,C3,D3] = bde3mep(b,b0,p3,q3,r3,s3,t3,u3,bc3,n3,opts);
