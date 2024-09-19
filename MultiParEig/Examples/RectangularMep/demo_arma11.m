%DEMO_ARMA11   demo for ARMA11
%
% We find critial points for the ARMA(1,1) model
%
% This is Example 14 from M.E.Hochstenbach, T.Kosir, B.Plestenjak: On the 
% solution of rectangular multiparameter eigenvalue problems, arXiv 2212.01867

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 6.11.2022

y = [2.4130 1.0033 1.2378 -0.72191 -0.81745 -2.2918 0.18213 0.073557 0.55248 2.0180 2.6593 1.1791];
y = y(:);

opts = [];
opts.showrank = 1;
% opts.heuristic = 0;

[points,val,err,cand] = arma11(y,opts);

N = length(y);
M = 800;
aset = linspace(-1,1,M);
gset = linspace(-1,1,M);
Z = zeros(M,M);
for i = 1:M
    for j = 1:M
        TC = diag(gset(j)*ones(N-1,1))+diag(ones(N-2,1),1); TC(N-1,N) = 1;
        TA = diag(aset(i)*ones(N-1,1))+diag(ones(N-2,1),1); TA(N-1,N) = 1;
        e = pinv(TC)*TA*y;
        Z(j,i) = norm(e)^2;
    end
end

contour(aset,gset,Z,100)
hold on
alpha = real(points(:,1));
gamma = real(points(:,2));
plot(alpha,gamma,'r*');
hold off
title('objective function \sigma^2','FontSize',14)
xlabel('\alpha_1','FontSize',14)
ylabel('\gamma_1','FontSize',14)

solution = [points val]

