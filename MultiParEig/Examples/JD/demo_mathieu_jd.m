%DEMO_MATHIEU_JD   demo for method twopareigs_jd
%
% This example computes first eigenmodes for the Mathieu two-parameter eigenvalue problem
% using Jacobi-Davidson method for two-parameter eigenvalue problems
%
% See also: TWOPAREIGS_JD, MATHIEU_MEP, TWOPAREIGS_IRA, TWOPAREIGS_SI

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 04.02.2018 Adapted to updated Jacobi-Davidson
% Last revision: 8.9.2015

n = 400;   % size of the matrices
neig = 20; % number of wanted eigenvalues
[A1,B1,C1,A2,B2,C2] = mathieu_mep(n,n,2,2,1); 

%% First version uses implicitly restarted Arnoldi with full vectors and
% Bartels-Stewart method for the Sylvester equation
tic; [lambda1,mu1,X1,Y1] = twopareigs_ira(A1,B1,C1,A2,B2,C2,neig); t1=toc
mu1

%% Second version uses subspace Arnoldi
opts = [];
opts.lowrank = neig;
opts.window = neig;
opts.arnsteps = 1;
opts.delta = eps*max(norm(A1),norm(A2));
opts.showinfo = 1;
opts.softlock = 0;
opts
tic; [lambda2,mu2,X1,Y1] = twopareigs_si(A1,B1,C1,A2,B2,C2,neig,opts); t2=toc
mu2

%% Third version is Jacobi-Davidson
% we multiply each equation by inverse of Ai
B1s = A1\B1; C1s = A1\C1; A1s = eye(n-2); 
B2s = A2\B2; C2s = A2\C2; A2s = eye(n-2);

opts = [];
opts.extraction = 'minmu';
opts.innersteps = -1;  % we solve correction equation exactly
opts.maxmaxsize = 20;
opts.maxsteps = 500;
opts.delta = 1e-10;
opts.forcereal = 1;
opts.refine = 3;
opts
tic; [lambda3,mu3,X3,Y3,X3l,Y3l,conv,hist] = twopareigs_jd(A1s,B1s,C1s,A2s,B2s,C2s,neig,opts); t3=toc
mu3
