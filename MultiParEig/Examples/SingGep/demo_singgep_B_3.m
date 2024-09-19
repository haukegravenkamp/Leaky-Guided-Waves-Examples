%DEMO_SINGGEP_B_3 Finite regular eigenvalues of a singular GEP are computed
% using a projection of normal rank
%
% This is Example 7.2 from 
% M.E. Hochstenbach, C. Mehl, B. Plestenjak: Solving singular generalized 
% eigenvalue problems part II: projection and augmentation, SIAM J. Matrix 
% Anal. Appl. 44 (2023) 1589-1618
% 
% For this singular pencil Matlab eig(A,B) return values that are all far 
% away from eigenvalue of (A,B), while singgep managed to compute 
% eigenvalues (projecttion method is used, but other two work well as well)
%
% See also SINGGEP

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Bor Plestenjak
% 27.06.2022

rng('default')

A = [1 -2 -100 0 0;1 0 -1 0 0; 0 0 0 1 -75; 0 0 0 0 2];
B = [zeros(4,1) eye(4)];
A(5,5) = 0;
B(5,5) = 0;

disp('Eigenvalues computed by eig(A,B) are completely wrong')
eigMatlab = eig(A,B)

opts = [];
opts.show = 1;

disp('Singgep (projection of normal rank) returns correct eigenvalues')
opts.method = 'project';
[lambda1,nrank1,Z1,d1,X1,Y1,U1,V1] = singgep(A,B,opts);
lambda1

disp('Singgep (rank-complete perturbation) returns correct eigenvalues')
opts.method = 'rank-complete';
[lambda2,nrank2,Z2,d2,X2,Y2,U2,V2] = singgep(A,B,opts);
lambda2

disp('Singgep (augmentation) returns correct eigenvalues')
opts.method = 'augment';
[lambda3,nrank3,Z3,d3,X3,Y3,U3,V3] = singgep(A,B,opts);
lambda3
