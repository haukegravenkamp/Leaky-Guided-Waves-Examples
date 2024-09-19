%DEMO_SINGGEP_RECT Finite regular eigenvalues of a singular rectangular pencil
% are computed using a projection of normal rank
% 
% This is Example 7.6 from 
% M.E. Hochstenbach, C. Mehl, B. Plestenjak: Solving singular generalized 
% eigenvalue problems part II: projection and augmentation, SIAM J. Matrix 
% Anal. Appl. 44 (2023) 1589-1618
%
% It originates from A. Emami-Naeini, P. M. Van Dooren: Computations of 
% zeros of linear multivariable systems, Automatica 18 (1982) 415-430.

A = [-2 -6 3 -7 6; 0 -5 4 -4 8; 0 2 0 2 -2; 0 6 -3 5 -6; 0 -2 2 -2 5];
B = [-2 7; -8 -5; -3 0; 1 5; -8 0];
C = [0 -1 2 -1 -1;1 1 1  0 -1; 0 3 -2 3 -1];
D = zeros(3,2);
E = eye(5);

F = [-A B;-C D];
G = [E zeros(5,2); zeros(3,7)];

opts = [];
opts.show = 1;
singgep(F,G,opts)



