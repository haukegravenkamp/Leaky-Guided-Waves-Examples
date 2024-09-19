%DEMO_SINGGEP_B_4 Values lambda such that A+lambda*B has a double 
% eigenvalue are computed as finite regular eigenvalues of a singular GEP 
% using a projection of normal rank
%
% This is Example 7.4 from 
% M.E. Hochstenbach, C. Mehl, B. Plestenjak: Solving singular generalized 
% eigenvalue problems part II: projection and augmentation, SIAM J. Matrix 
% Anal. Appl. 44 (2023) 1589-1618
%
% We take the double eigenvalue problems with two 20 x 20 matrices A,B and
% compute values lambda such that A+lambda*B has a double eigenvalue. They
% are eigenvalues of the singular pencil (Delta1,Delta0) of size 800 x 800.
% The method is compared to the method used in DEMO_SINGGEP4 that uses
% singular pencil of size 1200 x 1200.
%
% See also DEMO_SINGGEP4, SINGGEP

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Bor Plestenjak
% 27.06.2022

% We generate a random 20 x 20 example
N = 20;
rng('default')
A = randn(N);
B = randn(N);
I = eye(N);

% old approach from 2018 
disp('Approach that uses singular GEP of size 3n^2 x 3n^2')
[P,Q,R] = linearize_quadtwopar(A^2, A*B+B*A, 2*A, B^2, 2*B, I);
[Delta0,Delta1,Delta2] = twopar_delta(A,-B,-I,P,-Q,-R);
size_DeltaA = size(Delta0)
opts = [];
opts.method = 'rank-complete';
tic; [lambda1,nrank1,Z1,d1,X1,Y1,U1,V1,DA1,DB1] = singgep(Delta1,Delta0,opts); toc
found1 = length(lambda1)

% new approach from 2022
disp('New approach that uses singular GEP of size 2n^2 x 2n^2')
P = [A zeros(N); -I A];
Q = [B zeros(N); zeros(N) B];
R = [I zeros(N); zeros(N) I];
[Delta0,Delta1,Delta2] = twopar_delta(A,-B,-I,P,-Q,-R);
size_DeltaB = size(Delta0)
opts.method = 'project';
tic; [lambda2,nrank2,Z2,d2,X2,Y2,U2,V2] = singgep(Delta1,Delta0,opts); toc
found2 = length(lambda2)
