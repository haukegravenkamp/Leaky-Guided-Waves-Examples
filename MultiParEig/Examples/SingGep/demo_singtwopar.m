%DEMO_SINGTWOPAR Finite regular eigenvalues of a singular 2EP are computed
% using a rank-completing perturbation
%
% This is Example 7.1 from 
% M.E. Hochstenbach, C. Mehl, B. Plestenjak: Solving singular generalized 
% eigenvalue problems by a rank-completing perturbation, SIAM J. Matrix 
% Anal. Appl. 40 (2019) 1022-1046
%
% We solve system of bivariate polynomials
%
%  p1(x,y) = 1 + 2x + 3y + 4x^2 + 5xy + 6y^2 + 7x^3 + 8x^2y + 9xy^2 + 10y^3 = 0
%  p2(x,y) = 10 + 9x + 8y + 7x^2 + 6xy + 5y^2 + 4x^3 + 3x^2y + 2xy^2 +  y^3 = 0
% 
% using uniform determinantal representations with 5 x 5 matrices
% det(A1 + x*B1 + y*C1) = p1(x,y)
% det(A2 + x*B2 + y*C2) = p2(x,y)
% The corresponding two-parameter eigenvalue problems is singular and its
% finite regular eigenvalues are solutions (x,y) of p1(x,y)=0 and p2(x,y)=0
%
% See also DEMO_BIVARIATE

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Bor Plestenjak
% 23.05.2018

A1 = [ 0  0  4  1  0
       0  5  2  0  1
       6  3  1  0  0
       1  0  0  0  0
       0  1  0  0  0];
B1 = [ 0  0  7  0  0
       0  8  0 -1  0
       9  0  0  0 -1
       0  0  0  0  0
       0  0  0  0  0];
C1 = [ 0  0  0  0  0
       0  0  0  0  0
      10  0  0  0  0
       0 -1  0  0  0
       0  0 -1  0  0];   
A2 = [ 0  0  7  1  0
       0  6  9  0  1
       5  8 10  0  0
       1  0  0  0  0
       0  1  0  0  0];
B2 = [ 0  0  4  0  0
       0  3  0 -1  0
       2  0  0  0 -1
       0  0  0  0  0
       0  0  0  0  0];
C2 = [ 0  0  0  0  0
       0  0  0  0  0
       1  0  0  0  0
       0 -1  0  0  0
       0  0 -1  0  0];   
  
rng('default')

opts = [];
opts.method = 'rank-complete';
opts.show = 1;
opts.no_rotation = 1;
[x1,y1] = singtwopar(A1,-B1,-C1,A2,-B2,-C2,opts);
sol1 = [x1 y1]

% for comparison we solve this also with the staircase algorithm
opts = [];
opts.singular = 1;
[x2,y2] = twopareig(A1,-B1,-C1,A2,-B2,-C2,opts);
sol2 = [x2 y2]

% using projection to normal rank instead of rank completing perturbation
opts.method = 'project';
opts.show = 1;
opts.no_rotation = 1;
[x3,y3] = singtwopar(A1,-B1,-C1,A2,-B2,-C2,opts);
sol3 = [x3 y3]
