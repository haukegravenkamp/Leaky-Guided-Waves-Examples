%DEMO_SINGGEP_B_1 Finite regular eigenvalues of a singular GEP are computed
% using a projection of normal rank
%
% This covers Examples 6.2, 8.3 and 8.4 from 
% M.E. Hochstenbach, C. Mehl, B. Plestenjak: Solving singular generalized 
% eigenvalue problems part II: projection and augmentation, SIAM J. Matrix 
% Anal. Appl. 44 (2023) 1589-1618
% 
% 18 x 18 matrices A and B are constructed so that the KCF consists of
% blocks J_4(1), J_2(1), J_(1), N_2, N_1, L_2, L_1, L_2^2, and L_1^T
% normal rank of pencil (A,B) is 16
%
% Eigenvalues are computed three times using the projection of normal rank
% approach:
%   1) normal rank is computed accurately, we get correct eigenvalues
%   2) normal rank is underestimated, we get only three eigenvalue 1
%   3) normal rank is overestimated, we still get correct eigenvalues
%
% This example requires MCT Toolbox for Matlab

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Bor Plestenjak
% 27.06.2022

rng('default')

[A0,B0] = kcf(pstruct([2 1], [2 1], {[1,2,4]}, [1], [2,1]));
Q = randn(18);
Z = randn(18);
A = Q*A0*Z;
B = Q*B0*Z;

% we report 
origKCF = pguptri(A0,B0)

opts = [];
opts.show = 1;
opts.method = 'project';

% first run, normal rank is 16
disp('First run, normal rank is correctly estimated 16 in the algorithm')
[lambda1,nrank1,Z1,d1,X1,Y1,U1,V1] = singgep(A,B,opts);
nrank1
lambda1

% second run, normal rank is uderestimated to 15
disp('Second run, normal rank is underestimated to 15')
opts.nrank = 15;
[lambda2,nrank2,Z2,d2,X2,Y2,U2,V2] = singgep(A,B,opts);
lambda2

% first run, normal rank is 16
disp('Third run, normal rank is overestimated to 17')
opts.nrank = 17;
[lambda3,nrank3,Z3,d3,X3,Y3,U3,V3] = singgep(A,B,opts);
lambda3
