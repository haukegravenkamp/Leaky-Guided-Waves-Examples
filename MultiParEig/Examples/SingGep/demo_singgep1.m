%DEMO_SINGGEP1 Finite regular eigenvalues of a singular GEP are computed
% using a rank-completing perturbation
%
% This is Example 6.1 from 
% M.E. Hochstenbach, C. Mehl, B. Plestenjak: Solving singular generalized 
% eigenvalue problems by a rank-completing perturbation, SIAM J. Matrix 
% Anal. Appl. 40 (2019) 1022-1046
% 
% 7 x 7 matrices A and B are constructed so that the KCF consists of
% blocks J_1(1/2), J_1(1/3), N_1, L_1, and L_2^T
%
% See also DEMO_SINGGEP2, DEMO_SINGGEP3, DEMO_SINGGEP4, SINGGEP

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Bor Plestenjak
% 23.05.2018

A = [ -1    -1    -1    -1    -1    -1    -1
       1     0     0     0     0     0     0
       1     2     1     1     1     1     1
       1     2     3     3     3     3     3
       1     2     3     2     2     2     2
       1     2     3     4     3     3     3
       1     2     3     4     5     5     4];
   
B = [ -2    -2    -2    -2    -2    -2    -2
       2    -1    -1    -1    -1    -1    -1
       2     5     5     5     5     5     5
       2     5     5     4     4     4     4
       2     5     5     6     5     5     5
       2     5     5     6     7     7     7
       2     5     5     6     7     6     6];

opts = [];
opts.show = 1;
opts.method = 'rank-complete';

[lambda,nrank,Z,d,X,Y,U,V,DA,DB] = singgep(A,B,opts);

lambda

