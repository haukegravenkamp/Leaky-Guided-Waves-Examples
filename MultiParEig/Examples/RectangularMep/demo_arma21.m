%DEMO_ARMA21   demo for ARMA21
%
% We find critial points for the ARMA(2,1) model
%
% This is Example 16 from M.E.Hochstenbach, T.Kosir, B.Plestenjak: On the 
% solution of rectangular multiparameter eigenvalue problems, arXiv 2212.01867

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 6.11.2022

y = [0.41702 0.72032 0.01234 0.30233 0.14676 0.09234 0.18626];
y = y(:);

[points,val,err,cand] = arma21(y);

solution = [points val]
