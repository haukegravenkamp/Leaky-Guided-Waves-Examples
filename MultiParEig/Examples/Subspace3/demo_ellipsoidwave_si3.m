%DEMO_ELLIPSOIDWAVE_SI3 Ellipsoidal wave system is solved as a three-parameter 
% eigenvalue problem with Chebyshev collocation and subspace-iteration (example 3)
%
% This example returns Figure 5.2 (right) and third part of Table 5.2 in 
% M. Hochstenbach, K. Meerbergen, E. Mengi, B. Plestenjak: 
% Subspace methods for 3-parameter eigenvalue problems, arXiv 1802:07386
% 
% See also: DEMO_ELLIPSOIDWAVE_SI1, DEMO_ELLIPSOIDWAVE_SI2, DEMO_ELLIPSOIDWAVE_JD3

% We take the ellipsoidal wave problem discretized by Chebyshev collocation 
% using N = 300 points and we are looking for the eigenvalues with the 
% smallest |eta-1000| for the ellipsoid with semi-axes 1, 1.5, and 2 with a 
% fixed boundary (Dircichlet boundary condition) for the configuration 
% (rho,sigma,tau)=(0,0,0)

% We compute 80 eigenvalues using subspace iteration and the following settings:
%   - we shift the problem so that we are looking for smallest |eta|
%   - we limit delta matrices to size 15000 (we need even larger sizes as for target=200)
%   - we use 1 Arnoldi steps in the expansion
%   - we refine all Ritz values by one step of TRQI 
%   - we select Ritz values with the smallest |eta|
%   - we use "normalized" version of equations - where Ai = identity

% References: 
%   - M. Willatzen and L. C. Lew Yan Voon, Numerical implementation of the 
%     ellipsoidal wave equation and application to ellipsoidal quantum 
%     dots, Comput. Phys. Commun. 171 (2005) 1-18.
%   - B. Plestenjak, C. I. Gheorghiu, M. E. Hochstenbach: Spectral collocation 
%     for multiparameter eigenvalue problems arising from separable
%     boundary value problems, J. Comp. Phys. 298 (2015) 585-601
%   - M. Hochstenbach, K. Meerbergen, E. Mengi, B. Plestenjak: Subspace
%     methods for 3-parameter eigenvalue problems, arXiv 1802:07386

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Bor Plestenjak 26.02.2018

x0 = 1; y0 = 1.5; z0 = 2;
rho = 0; sigma = 0; tau = 0;
N = 300; % number of collocation nodes 
[A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,t1,t2,t3] = ellipsoidwave_mep(N,N,N,x0,y0,z0,rho,sigma,tau);

shifteta = 1000; % we are looking for eigenvalues such that [eta-shifteta| is minimal
shiftlambda = -5; 

neig = 80;

A1s = A1 - shiftlambda*B1 - shifteta*D1;
A2s = A2 - shiftlambda*B2 - shifteta*D2;
A3s = A3 - shiftlambda*B3 - shifteta*D3;

% we "normalize" matrices by multiplying all matricss with inv(A_i) 
n1 = size(A1s,1);
n2 = size(A2s,1);
n3 = size(A3s,1);
B1n = inv(A1s)*B1;
C1n = inv(A1s)*C1;
D1n = inv(A1s)*D1;
B2n = inv(A2s)*B2;
C2n = inv(A2s)*C2;
D2n = inv(A2s)*D2;
B3n = inv(A3s)*B3;
C3n = inv(A3s)*C3;
D3n = inv(A3s)*D3;
A1n = eye(n1);
A2n = eye(n2);
A3n = eye(n3);

opts = [];
opts.maxsteps = 10;
opts.arnsteps = 1; % for mildly interior eigenvalues we can not use zero Arnoldi steps anymore
opts.lowrank = 6; 
opts.window = 50;
opts.delta = 1e-8; 
opts.switcheps = 1e-2; 
opts.usesparse = 0;
opts.showinfo = 1;
opts.svdfilter = 1e-5;
opts.maxdetsize = 15000; 
opts.tol = eps;
opts.etamax = 25; % we limit search so that only candidates with |eta-1000|<25 are tested
opts.refine1 = 3;  % we also need more refinement steps for mildly interor eigenvalues  
opts.refine2 = 3;
opts

% we want the same example
if verLessThan('matlab', '7.7')
    rand('state',1);
else
    rng('default')
    rng(0);
end

tic
fprintf('Computing eigenvalues (lambda,mu,eta) for N=%d with smallest |eta-%d| using subspace iteration...\n',N,shifteta)
[lambda,mu,eta,X1,X2,X3,flag,hist] = threepareigs_si(A1n,B1n,C1n,D1n,A2n,B2n,C2n,D2n,A3n,B3n,C3n,D3n,neig,opts);
t1 = toc

lambda = lambda + shiftlambda;
eta = eta + shifteta;

plot(abs(eta-shifteta),'.','MarkerSize',20)

[tmp,ord] = sort(abs(eta-shifteta));

% solution (20 eigenvalues closest to the target)
ind = ord(1:20);
real([ind lambda(ind) mu(ind) eta(ind)])


