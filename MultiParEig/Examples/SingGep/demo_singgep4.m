%DEMO_SINGGEP4 Finite regular eigenvalues of a singular GEP are computed
% using a rank-completing perturbation
%
% This is Example 6.4 from 
% M.E. Hochstenbach, C. Mehl, B. Plestenjak: Solving singular generalized 
% eigenvalue problems by a rank-completing perturbation, SIAM J. Matrix 
% Anal. Appl. 40 (2019) 1022-1046
%
% We take the double eigenvalue problems with two 10 x 10 matrices A,B and
% compute values lambda such that A+lambda*B has a double eigenvalue. They
% are eigenvalues of the singular pencil (Delta1,Delta0) of size 300 x 300
% In this example the staircase algorithm in MultiParEig fails to extract
% finite regular eigenvalues
%
% The KCF contains 100 N1, 5 L_4^T, 5 L_5^T, 5 L_5, 5 L_6 and 90 J1 blocks.
%
% See also DEMO_DOUBLE_EIG_MP, SINGGEP, DEMO_SINGGEP1, DEMO_SINGGEP2, DEMO_SINGGEP3

% Reference: A. Muhic, B. Plestenjak:  A method for computing all values 
% lambda such that A + lambda*B has a multiple eigenvalue, Linear Algebra 
% Appl. 440 (2014) 345-359.

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Bor Plestenjak
% 23.05.2018

% We generate a random 10 x 10 example such that the staircase method fails
rng('default')
rng(3) 
A = rand(10);
B = rand(10);
I = eye(10);

[P,Q,R] = linearize_quadtwopar(A^2, A*B+B*A, 2*A, B^2, 2*B, I);
[Delta0,Delta1,Delta2] = twopar_delta(A,-B,-I,P,-Q,-R);

opts = [];
opts.singular = 1;

% staircase method fails to extract eigenvalues
lambda1 = twopareig(A,-B,-I,P,-Q,-R,opts)
found1 = size(lambda1,1)

% singgep finds all 90 finite eigenvalues

opts = [];
opts.show = 1;
opts.method = 'rank-complete';

[lambda2,nrank,Z,d,X,Y,U,V,DA,DB] = singgep(Delta1,Delta0,opts);
found2 = size(lambda2,1)

