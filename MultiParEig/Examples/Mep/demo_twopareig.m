%DEMO_TWOPAREIG   demo nonsingular two-parameter eigenvalue problems with 2 x 2 matrices
%
% We solve a two-parameter eigenvalue problem
%
% A1 x = lambda B1 x + mu C1 x 
% A2 y = lambda B2 y + mu C2 y,
%
% where 
%
% A1 = [1  2; 3  4]; B1 = [3  1; -1 1]; C1 = [2  1; 5 1],
% A2 = [1 -2; 3 -5]; B2 = [1 -1; -2 3]; C2 = [2 -1; 3 1].
%
% The output should include:
%
% eigenvalues =
%
%   -3.5718             5.6063          
%    3.9014            -1.0824          
%   -0.1364 + 0.0800i   0.0259 + 0.2820i
%   -0.1364 - 0.0800i   0.0259 - 0.2820i
%
% See also: TWOPAREIG

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

% matrices of the two-parameter eigenvalue problem
A1 = [1  2;  3  4]; 
B1 = [3  1; -1  1]; 
C1 = [2  1;  5  1];
A2 = [1 -2;  3 -5]; 
B2 = [1 -1; -2  3]; 
C2 = [2 -1;  3  1];

% Delta0 operator determinant
Delta0 = kron(B1,C2) - kron(C1,B2)

% rank of matrix Delta0 is 4 -> Delta0 is nonsingular
ran = rank(Delta0)

% we solve the two-parameter eigenvalue problem
[lambda,mu,Xr,Yr,Xl,Yl] = twopareig(A1,B1,C1,A2,B2,C2);

eigenvalues = [lambda mu]

% check that eigenvalues and eigenvectors are correct
for k = 1:size(eigenvalues,1)
    minsing1 = min(svd((A1-lambda(k)*B1-mu(k)*C1)));
    minsing2 = min(svd((A2-lambda(k)*B2-mu(k)*C2)));
    normres1 = norm((A1-lambda(k)*B1-mu(k)*C1)*Xr(:,k));
    normres2 = norm((A2-lambda(k)*B2-mu(k)*C2)*Yr(:,k));
    normres3 = norm(Xl(:,k)'*(A1-lambda(k)*B1-mu(k)*C1));
    normres4 = norm(Yl(:,k)'*(A2-lambda(k)*B2-mu(k)*C2));
    fprintf('Minimal singular values : (%7.1e, %7.1e), right residuals: (%7.1e, %7.1e), left residuals: (%7.1e,%7.1e)\n',...
        minsing1,minsing2,normres1,normres2,normres3,normres4)
end
    

