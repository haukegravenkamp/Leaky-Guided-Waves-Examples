function L = quad_left_compression_matrix(n)

%QUAD_LEFT_COMPRESSION_MATRIX Left compression matrix for quadratic rectangular 2EP
%
% L = quad_left_compression_matrix(n) returns left compression matrix for 
% the quadratic rectangular twoparameter eigenvalue problem 
%
% (A00 + l*A10 + u*A01 + l^2*A20 + l*u*A11 + u^2*A02)*x = 0 
%
% where matrices A_ij are of size (n+1) x n

% Bor Plestenjak, 2022

% See: M.E.Hochstenbach, T.Kosir, B.Plestenjak: On the solution of 
% rectangular multiparameter eigenvalue problems, arXiv 2212.01867

tmp = (1:3*n+1)';
e = ones(3*n+1,1);
P = [kron(tmp,e) kron(e,tmp)];
L = zeros(3*n*(n+1),(3*n+1)^2);
ind = 0;
for k = 1:(3*n+1)^2
    x = P(k,1);
    y = P(k,2);
    if x<=n+1 && y<=n+1 && x<y % situation (y_j,y_k), j<k
        ind = ind + 1;
        L(ind,k) = 1;
    elseif x<=n+1 && y>n+1 % situation (y_j,s_q) or (y_j,t_q)
        ind = ind + 1;
        L(ind,k) = 1;
    elseif n+1<x && x<2*n+2 && 2*n+1<y && x<=y-n % situation (s_p,t_q), p<q
        ind = ind + 1;
        L(ind,k) = 1;
    end
end
