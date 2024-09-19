function [A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,lambda,mu,eta] = random_sparse_3ep(n,dens)

%RANDOM_SPARSE_3EP   random sparse three-parameter eigenvalue problem
%
%[A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,lambda,mu,eta] = RANDOM_SPARSE_3EP(n,dens) builds
% a random three-parameter eigenvalue problem
%   A1 x = l B1 x + u C1 x + e D1 x
%   A2 y = l B2 y + u C2 y + e D2 y
%   A3 z = l B3 z + u C3 z + e D3 z
% with sparse matrices of size n x n and calculates eigenvalues (lambda,mu,eta)
%
% Input:
%   n : size of matrices
%   dens : density parameter used in the construction
%
% The obtained three-parameter eigenvalue problem is similar to a 
% three-parameter eigenvalue problem with diagonal matrices and this 
% enables us to compute all the eigenvalues and use them for tests
% in numerical methods

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% References: 
%  - M. Hochstenbach, K. Meerbergen, E. Mengi, B. Plestenjak: Subspace
%    methods for 3-parameter eigenvalue problems, arXiv 1802:07386

% Last revision 26.02.2018

% we want the same example
if verLessThan('matlab', '7.7')
    rand('state',1);
else
    rng('default')
    rng(0);
end

U1 = 0.3*sprand(n,n,dens)+speye(n);
U2 = 0.3*sprand(n,n,dens)+speye(n);
U3 = 0.3*sprand(n,n,dens)+speye(n); 
V1 = 0.3*sprand(n,n,dens)+speye(n);
V2 = 0.3*sprand(n,n,dens)+speye(n); 
V3 = 0.3*sprand(n,n,dens)+speye(n);
a1 = rand(n,1)-0.5;
b1 = rand(n,1)+2;
c1 = rand(n,1);
d1 = rand(n,1)-1;
a2 = rand(n,1)-0.5;
b2 = rand(n,1);
c2 = rand(n,1)+2;
d2 = rand(n,1)+0.5;
a3 = rand(n,1)-0.5;
b3 = rand(n,1)-1;
c3 = rand(n,1);
d3 = rand(n,1)+2;

lambda = zeros(n^3,1);
mu = zeros(n^3,1);
eta = zeros(n^3,1);

% in order to obtain all eigenvalues we need to solve two-parameter problem with diagonal matrices
% eigenvalues are intersections of lines a1+x b1+y c2=0 and a2+x b2+y c3=0

for k=1:n
   for j=1:n
       for r=1:n
            tmp = [b1(k),c1(k),d1(k); b2(j),c2(j),d2(j); b3(r),c3(r),d3(r)]\[a1(k);a2(j);a3(r)];
            indeks = (k-1)*n^2 + (j-1)*n + r;
            lambda(indeks) = tmp(1);
            mu(indeks) = tmp(2);
            eta(indeks) = tmp(3);
       end
   end
end

% we multiply matrices by nonsingular matrices so that they are not diagonal
A1 = V1*spdiag(a1)*U1;
B1 = V1*spdiag(b1)*U1;
C1 = V1*spdiag(c1)*U1;
D1 = V1*spdiag(d1)*U1;
A2 = V2*spdiag(a2)*U2;
B2 = V2*spdiag(b2)*U2;
C2 = V2*spdiag(c2)*U2;
D2 = V2*spdiag(d2)*U2;
A3 = V3*spdiag(a3)*U3;
B3 = V3*spdiag(b3)*U3;
C3 = V3*spdiag(c3)*U3;
D3 = V3*spdiag(d3)*U3;

end

function A = spdiag(d)

n = length(d);
A = spconvert([1:n; 1:n; d']');

end
