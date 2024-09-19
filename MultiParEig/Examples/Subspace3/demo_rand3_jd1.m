%DEMO_RAND3_JD1 Random three-parameter eigenvalue problem is solved
% with Jacobi-Davidson for target eta=-0.8 (example 1)
%
% See also: DEMO_RAND3_JD2, DEMO_RAND3_JD3
%
% This example computes 50 eigenvalues and plots Figure 5.5 in 
% M. Hochstenbach, K. Meerbergen, E. Mengi, B. Plestenjak: 
% Subspace methods for 3-parameter eigenvalue problems, arXiv 1802:07386
%
% We use a random generated 3-parameter eigenvalue problem with 100 x 100 
% matrices and compute eigenvalues with the smallest |eta+0.8|

% We compute 50 eigenvalues using Jacobi-Davidson and the following settings:
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

% we generate a random 3-parameter eigenvalue problem 
N = 100; % size of matrices
disp('Generating a random 3-parameter eigenvalue problem...')
[A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,lambda0,mu0,eta0] = random_sparse_3ep(N,0.04); 
A1 = full(A1); B1 = full(B1); C1 = full(C1); D1 = full(D1);
A2 = full(A2); B2 = full(B2); C2 = full(C2); D2 = full(D2);
A3 = full(A3); B3 = full(B3); C3 = full(C3); D3 = full(D3);

% Figure 
plot(lambda0,eta0,'.')

shifteta = -0.8; % we are looking for eigenvalues such that [eta-shift| is minimal

% we sort the eigenvalues by eta
[tmp,ord] = sort(abs(eta0-shifteta));
lambdas = lambda0(ord);
etas = eta0(ord);
mus = mu0(ord);

neig = 50;

A1s = A1 - shifteta*D1;
A2s = A2 - shifteta*D2;
A3s = A3 - shifteta*D3;

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
opts.maxsteps = 500;   % maximum number of outer iterations
opts.extraction = 'mineta'; 
opts.reschange = 0; % 10^(-6); % we change to minimal residual when residual is less then epschange (0 as we do not use this here)
opts.innersteps = -1;   % number of GMRES steps (-1 = exact solve correction equation)
opts.innertol = 1e-15;  % tolerance for the GMRES method
opts.target = [0 0 0];
opts.delta = 1e-10; 
opts.switcheps = 1e-6; 
opts.showinfo = 2; % set to 2 to see less information
opts.harmonic = 0;
opts.window = 0;
opts.refine = 3; % 
opts.forcereal = 1;
opts.selcrit1 = 1e-1;
opts.selcrit2 = 1e-4;
opts

tic
disp('Computing 50 eigenvalues (lambda,mu,eta) with smallest |eta+0.8| for the random 3-parameter problem with N=100 using J-D...')
[lambda,mu,eta,X1,X2,X3,Y1,Y2,Y3,res,hist] = threepareigs_jd(A1n,B1n,C1n,D1n,A2n,B2n,C2n,D2n,A3n,B3n,C3n,D3n,neig,opts);
toc

eta = eta + shifteta;

[tmp,ord] = sort(abs(eta-shifteta));

% now we compute indices of the computed eigenvalues
n1 = length(lambda);
indices = [];
for k=1:n1
    dist = abs(lambdas-lambda(k)).^2 + abs(mus-mu(k)).^2 + abs(etas-eta(k)).^2;
    [tmp,pos] = min(dist);
    indices(k) = pos;
end

% indices of the computed eigenvalues (how close are they to the target) in the order of retrieval
indices