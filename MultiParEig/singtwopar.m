 function [lambda,mu] = singtwopar(A1,B1,C1,A2,B2,C2,opts)

%SINGTWOPAR Solve a singular two-parameter eigenvalue problem using
% solver for singular generalized eigenvalue problems
% 
% [lambda,mu] = SINGTWOPAR(A1,B1,C1,A2,B2,C2) returns
% eigenvalues of a (singular) two-parameter eigenvalue problem
%
% A1 x = lambda B1 x + mu C1 x 
% A2 y = lambda B2 y + mu C2 y
%
% Input:
%   - A1, B1, C1, A2, B2, C2: matrices
%
% Output: 
%   - lambda, mu: eigenvalue parts (eigenvalues are (lambda(j),mu(j))
%
% Options in opts:
%   - sc_steps (0): number of initial staircase steps to reduce the size
%     of the problem. Reduction is stopped if Delta matrices become 
%     rectangular
%   - no_rotation (0): do not use random rotation of eigenvalues that
%     prevents that different eigenvalues have equal components
%   - method ('project'): which method from singgep should be used
%   - all options available for staircase_step_cr_np and singgep

% In 2022 the method has been updated and modified to be able to deal better 
% with multiple eigenvalues and different eigenvalues having the same lambda 
% or mu component

% Reference: Algorithm 2 in 
% M.E. Hochstenbach, C. Mehl, B. Plestenjak: Solving singular generalized 
% eigenvalue problems by a rank-completing perturbation, SIAM J. 
% Matrix Anal. Appl. 40 (2019) 1022-1046

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Bor Plestenjak
% 26.06.2022

if nargin<7, opts=[];  end

if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(A1,B1,C1,A2,B2,C2);
end

if ~isfield(opts,'method'),      opts.method = 'project';  end
if isfield(opts,'sc_steps'),     sc_steps = opts.sc_steps;       else, sc_steps=0;                   end
if isfield(opts, 'no_rotation'), no_rotation = opts.no_rotation; else, no_rotation = 0;              end

[Delta0,Delta1,Delta2] = twopar_delta(A1,B1,C1,A2,B2,C2);

% random rotation of parameters to ensure that different eigenvalues do not
% share the same lambda or mu component (unless eigenvalue is multiple)
if no_rotation
    fi = numeric_t(0,class_t);
else
    fi = rand(1,class_t)*numeric_t(pi,class_t);
end
c = cos(fi); s = sin(fi);
Delta1r =  c*Delta1 + s*Delta2;
Delta2r = -s*Delta1 + c*Delta2;
B1r = c*B1 + s*C1;
C1r = -s*B1 + c*C1;
B2r = c*B2 + s*C2;
C2r = -s*B2 + c*C2;

% optional initial staircase reduction of singular pencils 
% for some applications it can reduce the size of the problem, for other
% applications it might not work at all
if sc_steps>0
    OldDelta = {Delta0,Delta1r,Delta2r};
    for k=1:sc_steps
        Delta = staircase_step_cr_np(OldDelta,opts);
        [n1,n2] = size(Delta{1});
        if n1 ~= n2
            Delta = OldDelta;
            break
        else
            OldDelta = Delta;
        end
    end
    Delta0 = Delta{1};  
    Delta1r = Delta{2};  
    Delta2r = Delta{3};
end

% We solve a singular two-parameter eigenvalue problem by first applying
% singgep to the pencil (Delta1,Delta0) and then by inserting computed
% lambda's in the pencil of the first and the second equation and comparing 
% the computed mu's. We assume that to each lambda there corresponds exactly 
% one mu.
lambda = singgep(Delta1r,Delta0,opts);
n = length(lambda);
sol = [];
for k = 1:n
    M1 = A1 - lambda(k)*B1r;
    M2 = A2 - lambda(k)*B2r;
    mu1 = singgep(M1,C1r);
    mu2 = singgep(M2,C2r);
    mu_part = closest(mu1,mu2); % we compare mu1 and mu2 and take the closest pair
    sol = [sol; lambda(k) mu_part];
end

lambda = sol(:,1);
mu = sol(:,2);

% inverse rotation of parameters to transform solutions into correct
% solutions of the initial problem
lambdar = lambda; mur = mu;
lambda = c*lambdar - s*mur;
mu = s*lambdar + c*mur;
end

function x = closest(a,b)
    n = length(a);
    m = length(b);
    M = abs(a*ones(1,m)-ones(n,1)*(b.'));
    [tmp, indc] = min(M);
    [~, indr] = min(tmp);
    x = a(indc(indr));
end