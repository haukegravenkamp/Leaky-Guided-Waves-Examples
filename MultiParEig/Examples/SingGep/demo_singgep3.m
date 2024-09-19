%DEMO_SINGGEP3 Finite regular eigenvalues of a singular GEP are computed
% using a rank-completing perturbation
%
% This is Example 6.3 from 
% M.E. Hochstenbach, C. Mehl, B. Plestenjak: Solving singular generalized 
% eigenvalue problems by a rank-completing perturbation, SIAM J. Matrix 
% Anal. Appl. 40 (2019) 1022-1046
% 
% It originates from A. Edelman, Y. Ma: Staircase failures explained 
% by orthogonal versal forms, SIAM J. Matrix Anal. Appl. 21 (2000), 1004-1025.
% The 3 x 4 pencil has the structure J_2(0) and L_1
%
% See also DEMO_SINGGEP1, DEMO_SINGGEP2, DEMO_SINGGEP4, SINGGEP

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Bor Plestenjak
% 23.05.2018

d = 1.5e-8;
A = [0 0 1 0;0 0 0 1;0 0 0 0];
B = [d 0 0 0;0 d 0 0;0 0 1 0];

opts = [];
opts.show = 1;
opts.method = 'rank-complete';

% in the old version, the method failed to find 0 in unperturbed problems because eigenvalue 0 is not simple
% after the 2022 update, 0 is recognized as a finite eigenvalue 
[lambda,nrank,Z,d,X,Y,U,V,DA,DB] = singgep(A,B,opts);
lambda

% when we perturb the matrices using 1e-12 perturbation, guptri fails to find
% J_2(0), instead it returns J_1(0), for 1e-14 perturbation if does not return any regular eigenvalues 

[m,n] = size(A);
pert = 1e-14;
E = pert*rand(m,n);
F = pert*rand(m,n);
A1 = A + E;
B1 = B + F;

[lambda1,nrank1,Z1,d1,X1,Y1,U1,V1,DA1,DB1] = singgep(A1,B1,opts);
lambda1

pert = 1e-11;
E = pert*rand(m,n);
F = pert*rand(m,n);
A2 = A + E;
B2 = B + F;

[lambda2,nrank2,Z2,d2,X2,Y2,U2,V2,DA2,DB2] = singgep(A2,B2,opts);
lambda2