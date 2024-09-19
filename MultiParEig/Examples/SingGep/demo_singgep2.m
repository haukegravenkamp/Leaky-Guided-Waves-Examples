%DEMO_SINGGEP2 Finite regular eigenvalues of a singular GEP are computed
% using a rank-completing perturbation
%
% This is Example 6.2 from 
% M.E. Hochstenbach, C. Mehl, B. Plestenjak: Solving singular generalized 
% eigenvalue problems by a rank-completing perturbation, SIAM J. Matrix 
% Anal. Appl. 40 (2019) 1022-1046
% 
% It is based on example C3 from 
% J. Demmel, B. Kagström: Accurate solutions of ill-posed problems in 
% control theory, SIAM J. Matrix. Anal. Appl. 9 (1988), 126-145.
% The 4 x 5 pencil (A,B) has the KCF structure: L_2, J_1(1), J_1(2)
%
% See also DEMO_SINGGEP1, DEMO_SINGGEP3, DEMO_SINGGEP4, SINGGEP

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Bor Plestenjak
% 23.05.2018

A = [1 -2 -100 0 0;1 0 -1 0 0; 0 0 0 1 -75; 0 0 0 0 2];
B = [zeros(4,1) eye(4)];

opts = [];
opts.show = 1;
opts.method = 'rank-complete';

[lambda,nrank,Z,d,X,Y,U,V,DA,DB] = singgep(A,B,opts);
lambda

% we perturb the matrices using 1e-6 perturbation, so we have to use a
% different separation constant

[m,n] = size(A);
pert = 1e-6;
E = pert*rand(m,n);
F = pert*rand(m,n);
A1 = A + E;
B1 = B + F;
nrank = rank(rand*A1+rand*B1);

opts = [];
opts.show = 1;
opts.method = 'rank-complete';
opts.nrank = nrank;
opts.tau = 1e-2;
opts.sep = 1e-5;
opts.sepinf = 100*eps;

[lambda1,nrank1,Z1,d1,X1,Y1,U1,V1,DA1,DB1] = singgep(A1,B1,opts);
lambda1
