function [mu,X] = rect_quad_multipareig(A,opts)

%RECT_QUAD_MULTIPAREIG  Solve a rectangular quadratic multiparameter eigenvalue problem
%
% [mu,X] = RECT_QUAD_MULTIPAREIG(A,opts) returns eigenvalues and eigenvectors 
% of a rectangular quadratic k-parameter eigenvalue problem for k = 2,3,4
%
% For k=2 (A{1},...,A{6} are (n+1) x n matrices) the form is:
%
% (A{1} + mu(1)*A{2} + mu(2)*A{3} + mu(1)^2*A{4} + mu(1)*mu(2)*A{5} + mu(2)^2*A{6})x = 0
% 
% For k=3 (A{1},...A{10} are (n+2) x n matrices) the form is:
%
% (A{1} + mu(1)*A{2} + mu(2)*A{3} + mu(3)*A{4} + mu(1)^2*A{5} + mu(1)*mu(2)*A{6} 
%       + mu(1)*mu(3)*A{7} + mu(2)^2*A{8} + mu(2)*mu(3)*A{9} + mu(3)^2*A{10})x = 0
% 
% For k=4 (A{1},...A{15} are (n+3) x n matrices) the form is:
%
% (A{1} + mu(1)*A{2} + mu(2)*A{3} + mu(3)*A{4} + mu(4)*A{5} + mu(1)^2*A{6} + mu(1)*mu(2)*A{7} 
%       + mu(1)*mu(3)*A{8} + mu(1)*mu(4)*A{9} + mu(2)^2*A{10} + mu(2)*mu(3)*A{11} 
%       + mu(2)*mu(4)*A{12} + mu(3)^2*A{13} + mu(3)*mu(4)*A{14} + mu(4)^2*A{15})x = 0
%
% Input:
%   - A : cell of 6 (k=2), 10 (k=3) or 15 (k=4) matrices of size (n+k-1) x n
%   - opts : options
% 
% We use transformation to a MEP combined with Vandermonde compression
%
% Options in opts:
%   - all options for related method joint_delta_eig 
% Output:
%   - lambda : eigenvalues as rows
%   - X : matrix of right eigenvectors

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% P. Holoborodko, Advanpix LLC.
% FreeBSD License, see LICENSE.txt

% Bor Plestenjak 2023

% See: M.E.Hochstenbach, T.Kosir, B.Plestenjak: Numerical methods for rectangular 
% multiparameter eigenvalue problems, with applications to finding optimal 
% ARMA and LTI models. Numer Linear Algebra Appl. 2023; e2540

% Compression approach is based on: F. F. Alsubaie: H2 Optimal Model 
% Reduction for Linear Dynamic Systems and the Solution of Multiparameter 
% Matrix Pencil Problems, PhD, Imperial College London, 2019.

narginchk(1,2);

% Analyse user supplied options, if any
if nargin < 2, opts = []; end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(A{:});
end

nA = numel(A);
if nA == 6
    k = 2;
elseif nA == 10
    k = 3;
elseif nA == 15
    k = 4;
else
   error('A must be a cell of 6, 10 or 15 matrices of size (n+k-1) x k for k=2,3 or 4')
end

% Make sure all inputs are of the same numeric type.
for r=1:numel(A)
    if ~isa(A{r}, class_t)
         A{r} = numeric_t(A{r},class_t);
    end
end

[m,n] = size(A{1});
if m ~= n+k-1
   error('Matrices must be rectangular of size (n+k-1) x n')
end
for r = 2:numel(A)
    [m1,n1] = size(A{r});
    if any([m1 n1] ~= [m n])    
        error('Matrices must be rectangular of size (n+k-1) x n')
    end
end

% linearization to a linear rectangular MEP
ZX = zeros(n+k-1,n,class_t);
Z = zeros(n,class_t);
I = eye(n,class_t);
LA = cell(1,k+1);
if k==2
    LA{1} = -[A{1} A{2} A{3}; Z -I  Z; Z  Z -I];
    LA{2} =  [ZX   A{4} A{5}; I  Z  Z; Z  Z  Z];
    LA{3} =  [ZX   ZX   A{6}; Z  Z  Z; I  Z  Z];
elseif k==3
    LA{1} = -[A{1} A{2} A{3}  A{4}; Z -I  Z  Z; Z  Z -I  Z; Z  Z  Z -I];
    LA{2} =  [ZX   A{5} A{6}  A{7}; I  Z  Z  Z; Z  Z  Z  Z; Z  Z  Z  Z];
    LA{3} =  [ZX    ZX  A{8}  A{9}; Z  Z  Z  Z; I  Z  Z  Z; Z  Z  Z  Z];
    LA{4} =  [ZX    ZX  ZX   A{10}; Z  Z  Z  Z; Z  Z  Z  Z; I  Z  Z  Z];
else
    LA{1} = -[A{1} A{2}  A{3}   A{4}  A{5}; Z -I  Z  Z  Z; Z  Z -I  Z  Z; Z  Z  Z -I  Z; Z  Z  Z  Z -I];
    LA{2} =  [ZX   A{6}  A{7}   A{8}  A{9}; I  Z  Z  Z  Z; Z  Z  Z  Z  Z; Z  Z  Z  Z  Z; Z  Z  Z  Z  Z];
    LA{3} =  [ZX    ZX  A{10}  A{11} A{12}; Z  Z  Z  Z  Z; I  Z  Z  Z  Z; Z  Z  Z  Z  Z; Z  Z  Z  Z  Z];
    LA{4} =  [ZX    ZX     ZX  A{13} A{14}; Z  Z  Z  Z  Z; Z  Z  Z  Z  Z; I  Z  Z  Z  Z; Z  Z  Z  Z  Z];
    LA{5} =  [ZX    ZX     ZX     ZX A{15}; Z  Z  Z  Z  Z; Z  Z  Z  Z  Z; Z  Z  Z  Z  Z; I  Z  Z  Z  Z];
end

% Vandermonde compression for linearized quad RMEP
TK = speye((k+1)*n);
R = lambda_compress(k,k);
TR = right_compression_matrix(n,k);
pom = (1:n)*(k+1);
if k == 2 
    K = TK(:,[pom-2 pom-1 pom]);
    M = kron(kron(speye(3),K),speye(n));
elseif k == 3
    K = TK(:,[pom-3 pom-2 pom-1 pom]);
    K1 = kron(kron(speye(16),K),speye(n^2));
    K2 = kron(kron(speye(4),K),speye(4*n^2));
    K3 = kron(kron(speye(16*n),K),speye(n));
    M = K3*K2*K1;
else
    K = TK(:,[pom-4 pom-3 pom-2 pom-1 pom]);
    K1 = kron(kron(speye(125),K),speye(n^3));
    K2 = kron(kron(speye(25),K),speye(5*n^3));
    K3 = kron(kron(speye(5),K),speye(25*n^3));
    K4 = kron(kron(speye(125*n),K),speye(n^2));
    K5 = kron(kron(speye(25*n),K),speye(5*n^2));
    K6 = kron(kron(speye(125*n^2),K),speye(n));
    M = K6*K5*K4*K3*K2*K1;
end
    
PTR = M*kron(R,TR);
RWS = left_quad_compression(n,k);
NN = size(RWS,1);

DW = cell(k+1,1);
for j = 1:k+1
    DW{j} = zeros(NN,class_t);
end

SL = cell(k,k+1);
for i = 1:NN
    for j = 1:k
        for r = 1:k+1
            SL{j,r} = LA{r}(RWS(i,j),:);
        end
    end
    DeltaRow = multipar_delta(SL);
    for j = 1:k+1
        DW{j}(i,:) = DeltaRow{j}*PTR;
    end
end

opts.singular = 1;
mu = joint_delta_eig(DW,opts);
m = length(mu);

if nargout>1
    X = zeros(n,m);
    for j = 1:m
        if k==2
            W = A{1} + mu(j,1)*A{2} + mu(j,2)*A{3} + mu(j,1)^2*A{4}+ mu(j,1)*mu(j,2)*A{5} + mu(j,2)^2*A{6};
        elseif k==3
            W = A{1} + mu(j,1)*A{2} + mu(j,2)*A{3} + mu(j,3)*A{4} + mu(j,1)^2*A{5} + mu(j,1)*mu(j,2)*A{6} ... 
               + mu(j,1)*mu(j,3)*A{7} + mu(j,2)^2*A{8} + mu(j,2)*mu(j,3)*A{9} + mu(j,3)^2*A{10};
        else
            W = A{1} + mu(j,1)*A{2} + mu(j,2)*A{3} + mu(j,3)*A{4} + mu(j,4)*A{5} + mu(j,1)^2*A{6} + mu(j,1)*mu(j,2)*A{7} ...
               + mu(j,1)*mu(j,3)*A{8} + mu(j,1)*mu(j,4)*A{9} + mu(j,2)^2*A{10} + mu(j,2)*mu(j,3)*A{11} ...
               + mu(j,2)*mu(j,4)*A{12} + mu(j,3)^2*A{13} + mu(j,3)*mu(j,4)*A{14} + mu(j,4)^2*A{15};
        end
        X(:,j) = min_sing_vec(W,0);
    end
end

end

% ----------------------------------------------------------------------
% Auxiliary method for the compression of lambda part 
function R = lambda_compress(n,k)

tmp = [1 primes(100)];
x = tmp(1:n+1);
rows = x;

for j=2:k
    rows = kron(rows,x);
end

cols = unique(sort(rows));

R = spalloc(length(rows),length(cols),length(rows));
for i = 1:length(rows)
    [a,b] = find(cols == rows(i));
    R(i,b) = 1;
end

end

% ----------------------------------------------------------------------
% Auxiliary method for the compression from left size
function rows = left_quad_compression(n,k)

if k==2
    N = 3*n+1;
    nc = N^2;
    nr = nchoosek(n+1,2)*nchoosek(4,2);
    rowind = 0;
    rows = zeros(nr,2);
    
    % first group is (s1,s2), s1<s2
    for k1 = 1:n
        for k2 = k1+1:n+1
            col = kron_index([k1 k2],n);
            rowind = rowind +1; 
            rows(rowind,:)=[k1 k2];
        end
    end
    % second group is (s1,ap), (s1,bp)
    for k1 = 1:n+1
        for k2 = n+2:N
            col = kron_index([k1 k2],n);
            rowind = rowind +1; 
            rows(rowind,:)=[k1 k2];
        end
    end
    % third group is (ap,bq), p<=q
    for k1 = n+2:2*n+1
        for k2 = k1+n:N
            col = kron_index([k1 k2],n);
            rowind = rowind +1; 
            rows(rowind,:)=[k1 k2];
        end
    end
elseif k==3
    N = 4*n+2;
    nc = N^3;
    nr = nchoosek(n+2,3)*nchoosek(6,3);
    rowind = 0;
    rows = zeros(nr,3);
    
    % first group is (s1,s2,s3), s1<s2<s3
    for k1 = 1:n
        for k2 = k1+1:n+1
            for k3 = k2+1:n+2
                col = kron_index([k1 k2 k3],n);
                rowind = rowind +1; 
                rows(rowind,:)=[k1 k2 k3];
            end
        end
    end
    % second group is (s1,s2,ap), (s1,s2,bp), (s1,s2,cp), s1<s2
    for k1 = 1:n+1
        for k2 = k1+1:n+2
            for k3 = n+3:N
                col = kron_index([k1 k2 k3],n);
                rowind = rowind +1; 
                rows(rowind,:)=[k1 k2 k3];
            end
        end
    end
    % third group is (s1,ap,bq), (s1,ap,cq), (s1,bp,cq), p<=q
    for k1 = 1:n+2
        for k2 = n+3:3*n+2
            if k2<=2*n+2
                meja = [k2+n:3*n+2 k2+2*n:N];
            else
                meja = [k2+n:N];
            end
            for k3 = meja
                col = kron_index([k1 k2 k3],n);
                rowind = rowind +1; 
                rows(rowind,:)=[k1 k2 k3];
            end
        end
    end
    % last group is (ap,bq,cr), p<=q<=r
    for k1 = n+3:2*n+2
        for k2 = k1+n:3*n+2
            for k3 = k2+n:N
                col = kron_index([k1 k2 k3],n);
                rowind = rowind +1; 
                rows(rowind,:)=[k1 k2 k3];
            end
        end
    end    
else
    N = 5*n+3;
    nc = N^4;
    nr = nchoosek(n+3,4)*nchoosek(8,4);
    rowind = 0;
    rows = zeros(nr,4);
    
    % first group is (s1,s2,s3,s4), s1<s2<s3<s4
    for k1 = 1:n
        for k2 = k1+1:n+1
            for k3 = k2+1:n+2
                for k4 = k3+1:n+3
                    col = kron_index([k1 k2 k3 k4],n);
                    rowind = rowind +1; 
                    rows(rowind,:)=[k1 k2 k3 k4];
                end
            end
        end
    end
    % second group is (s1,s2,s3,ap), (s1,s2,s3,bp), (s1,s2,s3,cp), (s1,s2,s3,dp), s1<s2<s3
    for k1 = 1:n+1
        for k2 = k1+1:n+2
            for k3 = k2+1:n+3
                for k4 = n+4:N
                    col = kron_index([k1 k2 k3 k4],n);
                    rowind = rowind +1; 
                    rows(rowind,:)=[k1 k2 k3 k4];
                end
            end
        end
    end
    % third group is (s1,s2,ap,bq), (s1,s2,ap,cq), (s1,s2,ap,dq), (s1,s2,bp,cq), (s1,s2,bp,dq), (s1,s2,cp,dq), p<=q
    for k1 = 1:n+2
        for k2 = k1+1:n+3
            for k3 = n+4:4*n+3
                if k3<=2*n+3
                    meja4 = [k3+n:3*n+3 k3+2*n:4*n+3 k3+3*n:N]; % (ap,bq), (ap,cq), (ap, dq)
                elseif k3<=3*n+3
                    meja4 = [k3+n:4*n+3 k3+2*n:N]; % (bp,cq), (bp, dq)
                else
                    meja4 = [k3+n:N];
                end
                for k4 = meja4
                    col = kron_index([k1 k2 k3 k4],n);
                    rowind = rowind +1; 
                    rows(rowind,:)=[k1 k2 k3 k4];
                end
            end
        end
    end
    % fourth group is (s1,ap,bq,cr), (s1,ap,bq,dr), (s1,ap,cq,dr), (s1,bp,cq,dr), p<=q<=r
    for k1 = 1:n+3
        for k2 = n+4:3*n+3
            if k2<=2*n+3
                meja3 = [k2+n:3*n+3 k2+2*n:4*n+3]; % (ap,bq), (ap,cq)
            else
                meja3 = [k2+n:4*n+3];
            end
            for k3 = meja3
                if k3<=3*n+3
                    meja4 = [k3+n:4*n+3 k3+2*n:N];
                else
                    meja4 = [k3+n:N];
                end
                for k4 = meja4
                    col = kron_index([k1 k2 k3 k4],n);
                    rowind = rowind +1; 
                    rows(rowind,:)=[k1 k2 k3 k4];
                end
            end
        end
    end
    % fifth group is (ap,bq,cr,ds), p<=q<=r<=s
    for k1 = n+4:2*n+3
        for k2 = k1+n:3*n+3
            for k3 = k2+n:4*n+3
                for k4 = k3+n:N
                    col = kron_index([k1 k2 k3 k4],n);
                    rowind = rowind +1; 
                    rows(rowind,:)=[k1 k2 k3 k4];
                end
            end
        end
    end
end

end

% ----------------------------------------------------------------------
function ind = kron_index(row,n)

    k = length(row);
    N = (k+1)*n + k - 1;
    ind = row(1)-1;
    for j = 2:k
        ind = ind*N + row(j)-1;
    end
    ind = ind + 1;
end