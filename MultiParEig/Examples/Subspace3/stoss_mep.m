function [A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,t1,t2,t3] = stoss_mep(n1,n2,n3,opts)

%STOSS_MEP  Discretizes a 4-point boundary problem with three parameters as a three-parameter eigenvalue problem
%
% [A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,t1,t2,t3] = stoss_mep(n1,n2,n3)
%
% We have differential equation 
% 
% y''(x) + (lambda + 2*mu*cos(x) + 2*eta*cos(2x))*y(x)=0
%
% with boundary conditions y(0)=y(1)=y(2)=y(3)=0 
%
% We transform this into a 3-parameter eigenvalue problems
% 
% y1''(x1) + (lambda + 2 mu cos(x1) + 2 eta cos(2x1))y1(x1)=0 on [0,1] and b.c. y1(0)=y1(1)=0
% y2''(x2) + (lambda + 2 mu cos(x2) + 2 eta cos(2x2))y2(x2)=0 on [1,2] and b.c. y2(1)=y2(2)=0
% y3''(x3) + (lambda + 2 mu cos(x3) + 2 eta cos(2x3))y3(x3)=0 on [2,3] and b.c. y3(2)=y3(2)=0
%
% Input:
%   - n1, n2, n3: number of points for x1, x2, and x3
%   - opts : options
%
% Options in opts:
%   - fp_type: numeric type to use: 'single', 'double', or 'mp' (needs MCT),
%              default is the superior type of input data (x0,y0,z0)
% Output:
%   - A1,B1,C1,D1 : (n1-2) x (n1-2) matrices for the first equation
%   - A2,B2,C2,D2 : (n2-2) x (n2-2) matrices for the second equation
%   - A3,B3,C3,D3 : (n3-2) x (n3-2) matrices for the third equation
%   - t1, t2, t3 : ni points for 1st, 2nd, and 3rd equation (including endpoints)
%
% References: 
%  - H. J. Stoss, Ein Verfahren zur Berechnung des charakteristischen 
%    Exponenten der Differentialgleichung y''+(lambda + 2*lambda_1*cos(x)+2*lambda_2*cos(2x))y=0,
%    Numer. Math. 10 (1967) 423--436.
%  - M. Hochstenbach, K. Meerbergen, E. Mengi, B. Plestenjak: Subspace
%    methods for 3-parameter eigenvalue problems, arXiv 1802:07386

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Supports MCT

% Bor Plestenjak 26.02.2018

narginchk(3, 4);

if nargin<4; 
    opts = [];
end

if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = 'double';
end
opts.fp_type = class_t;

p = numeric_t('1',class_t);
q = numeric_t('0',class_t);
r = numeric_t('0',class_t);
s = numeric_t('-1',class_t);
t = @(x) -2*cos(x);
u = @(x) -2*cos(2*x);
bc = numeric_t('[1 0; 1 0]',class_t); 

% Eq. 1  (0 < t1 < 1)
[t1,A1,B1,C1,D1] = bde3mep(0,1,p,q,r,s,t,u,bc,n1,opts);
% Eq. 2  (1 < t2 < 2)
[t2,A2,B2,C2,D2] = bde3mep(1,2,p,q,r,s,t,u,bc,n2,opts);
% Eq. 3  (2 < t3 < 3)
[t3,A3,B3,C3,D3] = bde3mep(2,3,p,q,r,s,t,u,bc,n3,opts);


