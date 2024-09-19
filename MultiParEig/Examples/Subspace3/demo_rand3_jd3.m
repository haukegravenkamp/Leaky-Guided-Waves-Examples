%DEMO_RAND3_JD3 Random three-parameter eigenvalue problem is solved
% with Jacobi-Davidson for target (0,0,0) (example 2)
%
% See also: DEMO_RAND3_JD1, DEMO_RAND3_JD2
%
% We use a random generated 3-parameter eigenvalue problem with 100 x 100 
% matrices and compute eigenvalues with the smallest |eta|^2+|lambda|^2+|mu|^2. 

% We compute 50 eigenvalues using Jacobi-Davidson and the following settings:
%    - we use 5 GMRES steps to solve the correction equation aproximately
%    - as preonditioners we use inv(Ai-lambdaT*Bi-muT*Ci-etaT*Di), where (lambdaT,muT,etaT)=(0,0,0) is the target
%    - we use TRQI to compute an eigenvalue when Ritz pair has a small residual
%    - we select only Ritz values that satisfy the selection criteria

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

% we sort the eigenvalues by their distance from (0,0,0)
dist = abs(lambda0).^2 + abs(mu0).^2 + abs(eta0).^2;
[tmp,ord] = sort(dist);
lambdas = lambda0(ord);
etas = eta0(ord);
mus = mu0(ord);

neig = 50;
% parameters for the Jacobi-Davidson
opts = [];							
opts.minsize = 5;     
opts.maxsize = 10;
opts.maxsteps = 500;   % maximum number of outer iterations
opts.extraction = 'mindist'; 
opts.reschange = 0; % 10^(-6); % we change to minimal residual when residual is less then epschange (0 as we do not use this here)
opts.innersteps = 10; % -1;   % number of GMRES steps (-1 = exact solve correction equation)
opts.innertol = 1e-10;  % tolerance for the GMRES method
opts.target = [0 0 0];
opts.M1 = inv(A1);
opts.M2 = inv(A2);
opts.M3 = inv(A3);
opts.delta = 1e-10; 
opts.switcheps = 1e-6; % usually we use 1e4*opts.delta, this is a bit more agressive
opts.showinfo = 2; % set to 2 to see less information
opts.harmonic = 0;
opts.window = 0;
opts.refine = 3;
opts.forcereal = 1;
opts.selcrit1 = 1e-1;
opts.selcrit2 = 1e-4;
opts

tic
disp('Computing 50 eigenvalues (lambda,mu,eta) closest to (0,0,0) for random 3EP with N=100 using JD iteration...')
[lambda,mu,eta,X1,X2,X3,Y1,Y2,Y3,res,hist] = threepareigs_jd(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,neig,opts);
toc

% now we compute indices of the computed eigenvalues
n1 = length(lambda);
indices = [];
for k=1:n1
    dist = abs(lambdas-lambda(k)).^2 + abs(mus-mu(k)).^2 + abs(etas-eta(k)).^2;
    [tmp,pos] = min(dist);
    indices(k) = pos;
end

% indices of the computed eigenvalues
indices

