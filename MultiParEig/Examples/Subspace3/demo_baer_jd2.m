%DEMO_BAER_JD2 Baer wave system is solved as a three-parameter eigenvalue problem 
% with Chebyshev collocation and Jacobi-Davidson (example 2)
%
% This example returns last part of Table 5.5 in 
% M. Hochstenbach, K. Meerbergen, E. Mengi, B. Plestenjak: 
% Subspace methods for 3-parameter eigenvalue problems, arXiv 1802:07386
%
% See also: DEMO_BAER_JD1, DEMO_BAER_SI1, DEMO_BAER_SI2
%
% We are looking for wavenumbers closest to 10 of a body that is an 
% intersection of two elliptical paraboloids (we use paraboloidal 
% coordinates) with a fixed boundary (Dircichlet boundary condition)

% We discretize the Baer wave equations by Chebyshev collocation on 
% N = 300 points and we are looking for the eigenvalues with smallest 
% |eta-100|. Eigenfrequencies are square roots of parts eta of 
% eigenvalues (lambda,mu,eta).
%
% We compute 80 eigenvalues using Jacobi-Davidson and the following settings:
%    - we shift the problem so that we are looking for smallest |eta|
%    - we solve the correction equation exactly
%    - we use TRQI to compute an eigenvalue when Ritz pair has a small residual
%    - we select Ritz values with the smallest |eta|
%    - we select only Ritz values that satisfy the selection criteria
%    - we use "normalized" version of equations - where Ai = identity

% References: 
%  - M. Hochstenbach, K. Meerbergen, E. Mengi, B. Plestenjak: Subspace
%    methods for 3-parameter eigenvalue problems, arXiv 1802:07386

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Bor Plestenjak 26.02.2018

N = 300;
[A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,t1,t2,t3] = baer_mep(N,N,N,1,3,0,5);

shiftmu = 5;     % we have to shift Ai with Ci to make them nonsingular
shifteta = 100;  % eta = omega^2, we are looking for |eta| closest to 100

A1s = A1 - shiftmu*C1 - shifteta*D1;
A2s = A2 - shiftmu*C2 - shifteta*D2;
A3s = A3 - shiftmu*C3 - shifteta*D3;

% we normalize the matrices (multiply the equation with Ai inverse)
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

neig = 80;
opts = [];
opts.target = [0 0 0];
opts.M1 = [];
opts.M2 = [];
opts.M3 = [];
opts.harmonic = 0;
opts.extraction = 'mineta';
opts.innersteps = -1;  % we solve the correction equation exactly
opts.minsize = 5;
opts.maxsize = 10;
opts.maxsteps = 1000;
opts.delta = 1e-8;
opts.switcheps = 1e-1;
opts.window = 0;
opts.showinfo = 2;
opts.forcereal = 1;
opts.refine = 4;
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
[lambda,mu,eta,X1,X2,X3,Y1,Y2,Y3,convlog,hist] = threepareigs_jd(A1n,B1n,C1n,D1n,A2n,B2n,C2n,D2n,A3n,B3n,C3n,D3n,neig,opts);
toc

mu = mu + shiftmu;
eta = eta + shifteta;

Z1 = [X1; zeros(1,size(X1,2))];
Z2 = X2;
Z3 = [zeros(1,size(X3,2)); X3];

% we compute number of zeros on each interval as eigenvalues have unique indices (i,j,k)
cs = [];
for k = 1:length(lambda);
    cs(k,1) = count_sign_changes(Z1(:,k));
    cs(k,2) = count_sign_changes(Z2(:,k));
    cs(k,3) = count_sign_changes(Z3(:,k));
end

k = length(lambda);

[tmp,ord] = sort(abs(eta-shifteta));

% solution (20 eigenvalues closest to the target)
ind = ord(1:20);
real([ind lambda(ind) mu(ind) eta(ind) sqrt(eta(ind)) cs(ind,:)])

XX = [];
str = [];
tt = [t3; t2; t1];
for j = 1:6
  k = ind(j);  
  XX(:,j) = [Z3(:,k)*Z2(1,k)/Z3(end,k)*Z1(1,k)/Z2(end,k); Z2(:,k)*Z1(1,k)/Z2(end,k); Z1(:,k)];
  XX(:,j) = XX(:,j)/max(abs(XX(:,j)));
  str = [str; sprintf('\\omega_%d = %11.8f',j,sqrt(eta(k)))];
end
plot(tt,XX,'LineWidth',2)
title('Eigenfunctions for eigenfrequencies closest to 10')
axis([0 5 -1.1 1.1])
legend(cellstr(str),'Location','NorthEastOutside','FontSize',15)


