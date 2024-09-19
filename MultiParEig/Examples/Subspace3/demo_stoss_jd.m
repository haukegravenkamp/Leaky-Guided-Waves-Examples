%DEMO_STOSS_JD 4-point boundary problem is solved as a three-parameter eigenvalue problem 
% with Chebyshev collocation and Jacobi-Davidson 
%
% This example returns Figure 5.4 and Table 5.6 in 
% M. Hochstenbach, K. Meerbergen, E. Mengi, B. Plestenjak: 
% Subspace methods for 3-parameter eigenvalue problems, arXiv 1802:07386
%
% See also: STOSS_MEP

% We are looking for solutions (lambda,mu,eta) for Stoss differential equation 
% 
% y''(x) + (lambda + 2*mu*cos(x) + 2*eta*cos(2x))*y(x)=0
%
% with boundary conditions y(0)=y(1)=y(2)=y(3)=0 

% We discretize the equations by Chebyshev collocation on N = 200 points and 
% use Jacobi-Davidson to find eigenvales closest to (0,0,0) as these are
% the eigenvalues with the smallest indices.
%
% We compute 20 eigenvalues using Jacobi-Davidson and the following settings:
%    - we use 10 step of GMRES to solve the correction equation approximately
%    - we use inv(Ai) as preconditioners 
%    - we use TRQI to compute an eigenvalue when Ritz pair has a small residual
%    - we select Ritz values with the smallest |lambda|^2+|mu|^2+|eta|^2
%    - we select only Ritz values that satisfy the selection criteria

% References: 
%  - M. Hochstenbach, K. Meerbergen, E. Mengi, B. Plestenjak: Subspace
%    methods for 3-parameter eigenvalue problems, arXiv 1802:07386

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Bor Plestenjak 26.02.2018

%% First option - we use Chebyshev collocation on 12 points for comparison
%
% This gives the problems that is small enough that we can compute all
% eigenvalues. We return 20 eigenvalues closest to (0,0,0). For each
% eigenvalue we compute the corresponding index

disp('Computing eigenvalues (lambda,mu,eta) for N=12 close to (0,0,0) using full Delta matrices...')
[A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,t1,t2,t3] = stoss_mep(12,12,12);
[lambda1,mu1,eta1,W1,W2,W3] = threepareig(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3);

X1 = [zeros(1,size(W1,2)); W1; zeros(1,size(W1,2))];
X2 = [zeros(1,size(W2,2)); W2; zeros(1,size(W2,2))];
X3 = [zeros(1,size(W3,2)); W3; zeros(1,size(W3,2))];

cs1 = [];
for k=1:length(lambda1);
    cs1(k,1) = count_sign_changes(X1(:,k));
    cs1(k,2) = count_sign_changes(X2(:,k));
    cs1(k,3) = count_sign_changes(X3(:,k));
end

[tmp1,ord1] = sort(abs(lambda1).^2+abs(mu1).^2+abs(eta1).^2);

k = ord1(1:20); [k lambda1(k) mu1(k) eta1(k) cs1(k,:)]

%% We use Chebyshev collocation on 200 points to get more accurate results
%
% We use Jacobi-Davidson, using (0,0,0) as a target, preconditioning with
% inv(A1), inv(A2), inv(A3) and solve correction equations approximately
% using 10 steps of GMRES.
N = 200;
[A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,t1,t2,t3] = stoss_mep(N,N,N);

neig = 20;
opts = [];
opts.target = [0 0 0];
opts.M1 = inv(A1);
opts.M2 = inv(A2);
opts.M3 = inv(A3);
opts.harmonic = 0;
opts.extraction = 'mindist';
opts.innersteps = 10;
opts.minsize = 5;
opts.maxsize = 10;
opts.maxsteps = 200;
opts.delta = 1e-10;
opts.switcheps = 1e-6;
opts.window = 0;
opts.showinfo = 2;
opts.forcereal = 1;
opts.refine = 3;
opts

if verLessThan('matlab', '7.7')
    rand('state',1);
else
    rng('default')
    rng(0);
end

tic
fprintf('Computing eigenvalues (lambda,mu,eta) for N=%3d close to (0,0,0) using J-D...\n',N)
[lambda,mu,eta,X1,X2,X3,Y1,Y2,Y3,convlog,hist] = threepareigs_jd(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,neig,opts);
toc

Z1 = [zeros(1,size(X1,2)); X1; zeros(1,size(X1,2))];
Z2 = [zeros(1,size(X2,2)); X2; zeros(1,size(X2,2))];
Z3 = [zeros(1,size(X3,2)); X3; zeros(1,size(X3,2))];

cs = [];
for k=1:length(lambda);
    cs(k,1) = count_sign_changes(Z1(:,k));
    cs(k,2) = count_sign_changes(Z2(:,k));
    cs(k,3) = count_sign_changes(Z3(:,k));
end

k = length(lambda);

[(1:k)' lambda mu eta cs]

t1 = t1(end:-1:1);
t2 = t2(end:-1:1);
t3 = t3(end:-1:1);
Z1 = Z1(end:-1:1,:);
Z2 = Z2(end:-1:1,:);
Z3 = Z3(end:-1:1,:);

tt = [t1;t2;t3];
XX = [];
for col = 1:3
    for row = 1:3
        k = (row-1)*3+col;
        odv1a = (Z1(end,k)-Z1(end-1,k))/(t1(end)-t1(end-1));
        odv1b = (Z2(2,k)-Z2(1,k))/(t2(2)-t2(1));
        odv2a = (Z2(end,k)-Z2(end-1,k))/(t2(end)-t2(end-1));
        odv2b = (Z3(2,k)-Z3(1,k))/(t3(2)-t3(1));
        XX(:,k) = [Z1(:,k); Z2(:,k)*odv1a/odv1b; Z3(:,k)*odv1a/odv1b*odv2a/odv2b];
        XX(:,k) = XX(:,k) / max(abs(XX(:,k)));
        subplot(3,3,k)
        plot(tt,XX(:,k),'LineWidth',1)
        axis([0 3 -1.2 1.2])
        set(gca,'FontSize',10)
        title(sprintf('(%1d,%1d,%1d)',cs(k,1),cs(k,2),cs(k,3)),'FontSize',12)
    end
end

