%DEMO_SINGGEP_B_2 Finite regular eigenvalues of a singular GEP are computed
% using a) a projecton of normal rank, b) augmentation
%
% This is Example 7.1 from 
% M.E. Hochstenbach, C. Mehl, B. Plestenjak: Solving singular generalized 
% eigenvalue problems part II: projection and augmentation, SIAM J. Matrix 
% Anal. Appl. 44 (2023) 1589-1618
% 
% 7 x 7 matrices A and B are constructed so that the KCF consists of
% blocks J_1(1/2), J_1(1/3), N_1, L_1, and L_2^T
%
% See also DEMO_SINGGEP1, SINGGEP

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Bor Plestenjak
% 27.06.2022

rng('default')

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

disp('First method - projection of normal rank')
opts.method = 'project';
[lambda1,nrank1,Z1,d1,X1,Y1,U1,V1] = singgep(A,B,opts);
lambda1

disp('Second method - augmentation')
opts.method = 'augment';
[lambda2,nrank2,Z2,d2,X2,Y2,U2,V2] = singgep(A,B,opts);
lambda2

disp('Third method - augment A with U and V, and B with zero blocks')
opts.method = 'augment';
opts.DA=1;
opts.DB=0;
opts.SA=1;
opts.SB=0;
[lambda3,nrank3,Z3,d3,X3,Y3,U3,V3] = singgep(A,B,opts);
lambda3



