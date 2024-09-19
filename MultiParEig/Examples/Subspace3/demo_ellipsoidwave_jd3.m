%DEMO_ELLIPSOIDWAVE_JD3 Ellipsoidal wave system is solved as a three-parameter 
% eigenvalue problem with Chebyshev collocation and Jacobi-Davidson (example 3)
%
% This example returns Figure 5.1 (right) and third part of Table 5.2 in 
% M. Hochstenbach, K. Meerbergen, E. Mengi, B. Plestenjak: 
% Subspace methods for 3-parameter eigenvalue problems, arXiv 1802:07386
% 
% See also: DEMO_ELLIPSOIDWAVE_JD1, DEMO_ELLIPSOIDWAVE_JD2, DEMO_ELLIPSOIDWAVE_SI3

% We take the ellipsoidal wave problem discretized by Chebyshev collocation 
% using N = 300 points and we are looking for eigenvalues such that eta is 
% closest to 1000 for the ellipsoid with semi-axes 1, 1.5, and 2 with a 
% fixed boundary (Dircichlet boundary condition) for the configuration 
% (rho,sigma,tau)=(0,0,0).

% We compute 80 eigenvalues using Jacobi-Davidson and the following settings:
%    - we shift the problem so that we are looking for smallest |eta|
%    - we solve the correction equation exactly
%    - we use TRQI to compute an eigenvalue when Ritz pair has a small residual
%    - we select Ritz values with the smallest |eta|
%    - we select only Ritz values that satisfy the selection criteria
%    - we use "normalized" version of equations - where Ai = identity

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

shifteta = 1000; % we are looking for eigenvalues such that [eta-shift| is minimal
shiftlambda = -5; % we use this shift to make matrices Ai nonsingular

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

% parameters for the Jacobi-Davidson
opts = [];							
opts.minsize = 5; 
opts.maxsize = 10; 
opts.maxsteps = 1000; % maximum number of outer iterations
opts.extraction = 'mineta'; 
opts.reschange = 0; % we change to minimal residual when residual is less then epschange (o as we do not use this option)
opts.innersteps = -1;   % number of GMRES steps (-1 = solve correction equation exactly)
opts.innertol = 1e-15;  % tolerance for the GMRES method
opts.target = [0 0 0];
opts.delta = 1e-8; 
opts.switcheps = 1e-1; % usually we use 1e4*opts.delta, this is a bit more agressive 
opts.showinfo = 2; % set to 1 to see more information
opts.harmonic = 0;
opts.window = 0;
opts.refine = 4;
opts.forcereal = 1;
opts.selcrit1 = 1e-1;
opts.selcrit2 = 1e-4;
opts

% we want the same example
if verLessThan('matlab', '7.7')
    rand('state',1);
else
    rng('default')
    rng(0);
end

tic
disp('Computing eigenvalues (lambda,mu,eta) for N=300 with smallest |eta-1000| using JD iteration...')
[lambda,mu,eta,X1,X2,X3,Y1,Y2,Y3,res,hist] = threepareigs_jd(A1n,B1n,C1n,D1n,A2n,B2n,C2n,D2n,A3n,B3n,C3n,D3n,neig,opts);
toc

lambda = lambda + shiftlambda;
eta = eta + shifteta;

plot(abs(eta-shifteta),'.','MarkerSize',20)

[tmp,ord] = sort(abs(eta-shifteta));

% solution (20 eigenvalues closest to the target)
ind = ord(1:20);
real([ind lambda(ind) mu(ind) eta(ind)])

