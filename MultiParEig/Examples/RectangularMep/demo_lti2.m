%DEMO_LTI2   demo for LTI2
%
% We find critial points for the LTI(2) model
%
% This is Example 17 from M.E.Hochstenbach, T.Kosir, B.Plestenjak: On the 
% solution of rectangular multiparameter eigenvalue problems, arXiv 2212.01867

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 5.12.2022

y = [0.69582 0.68195 -0.24647 0.50437 -0.23207 0.34559 -0.19628 0.20553 -0.17737 0.11543];
y = y(:);

[points,val,err,cand] = lti2(y);

solution = [points val]
