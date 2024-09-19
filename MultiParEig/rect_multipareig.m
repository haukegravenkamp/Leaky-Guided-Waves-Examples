function [lambda,X] = rect_multipareig(A,opts)

%RECT_MULTIPAREIG   Solve a linear rectangular multiparameter eigenvalue problem
%
% [lambda,X] = rect_multipareig(A) returns eigenvalues and eigenvectors
% of a rectangular multiparameter eigenvalue problem
%
% A{1} x + lambda(1) A{2} x + ... + lambda(k) A{k+1} x = 0 
% 
% Input:
%   - A : cell array of size k+1 of matrices A{i}, all matrices have to be 
%         rectangular matrices of the same size (n+k-1) x n
%   - opts : options 
%
% Options in opts:
%   - method: 'compress' (default) or 'mep': method to use - 'compress' 
%     transforms problem into a joined system of GEPS, 'mep' transforms
%     problem into a k-parameter eigenvalue problem using k random
%     prejections
%   - fp_type: numeric type to use ('single', 'double', or 'mp' (needs MCT),
%     use only if you want to change the default - the superior type of input data
%   - delta (sqrt(eps)): treshold for eigenvalues of the rectangular
%     problem for the 'mep' method
%   - all options for multipareig (for 'mep' method) and joint_delta_eig
%     (for 'compress' method)
%
% Output:
%   - lambda : matrix m x k, each row is an eigenvalue
%   - X : matrix n x m with right eigenvectors

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% P. Holoborodko, Advanpix LLC.
% FreeBSD License, see LICENSE.txt

% Bor Plestenjak 2022
% BP 08.11.2023: Faster compress method using less memory

% See: M.E.Hochstenbach, T.Kosir, B.Plestenjak: Numerical methods for rectangular 
% multiparameter eigenvalue problems, with applications to finding optimal 
% ARMA and LTI models. Numer Linear Algebra Appl. 2023; e2540

% Compression method is based on: F. F. Alsubaie: H2 Optimal Model 
% Reduction for Linear Dynamic Systems and the Solution of Multiparameter 
% Matrix Pencil Problems, PhD, Imperial College London, 2019.

narginchk(1, 2);

% Parse user supplied options, if any
if nargin < 2, opts = []; end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(A{1});
end
if isfield(opts,'method'),    method    = opts.method;    else,   method = 'compress';                         end
if isfield(opts,'delta'),     delta     = opts.delta;     else,   delta = sqrt(eps(class_t));                  end

k = length(A) - 1;
for j = 1:k+1
    if ~isa(A{j},class_t), A{j} = numeric_t(A{j},class_t); end
end

[m,n] = size(A{1});
if m ~= n+k-1
   error('Matrices must be rectangular of size (n+k-1) x n')
end
for j = 2:k+1
   [m1,n1] = size(A{j});
   if [m1 n1] ~= [m n]
      error('Matrices must be rectangular of size (n+k-1) x n')
   end
end

if strcmp(method,'compress')
    TR = right_compression_matrix(n,k);
    [~,RWS] = left_compression_matrix(m,k);
    NN = size(RWS,1);

    DW = cell(1,k+1);
    for j = 1:k+1
        DW{j} = zeros(NN,class_t);
    end
        
    SL = cell(k,k+1);
    % we compute matrices DW row-by-row only the rows from RWS
    for r = 1:NN
        for i = 1:k
            SL{i,1} = A{1}(RWS(r,i),:);
            for j = 2:k+1
                SL{i,j} = -A{j}(RWS(r,i),:);
            end
        end
        DeltaRow = multipar_delta(SL);
        for j = 1:k+1
            DW{j}(r,:) = DeltaRow{j}*TR;
        end
    end

    lambda = joint_delta_eig(DW,opts);
    if nargout>1
        X = zeros(n,length(lambda),class_t);
        for i = 1:length(lambda)
            TMP = A{1};
            for j = 1:k
                TMP = TMP + lambda(i,j)*A{j+1};
            end
            X(:,i) = min_sing_vec(TMP,0); % we have to use SVD as matrices are rectangular
        end
    end
   
elseif strcmp(method,'mep')
    
    % we project original problem into a multparameter eigenvalue problem
    W = cell(k,k+1);
    for i = 1:k
        P = orth(randn(m,n,class_t)).'; % random projection matrix
        W{i,1} = P*A{1};
        for j = 2:k+1
            W{i,j} = -P*A{j};
        end
    end

    vecnorm = zeros(k+1,1);
    for j=1:k+1
        vecnorm(j) = norm(A{j},'fro');
    end

    % we solve the multiparameter eigenvalue problem that has extra solutions
    % and then return only solutions that solve the rectangular problem
    lambdaC = multipareig(W,opts);

    lambda = [];
    X = [];
    for j = 1:size(lambdaC,1)
        mat = A{1};
        ocena = vecnorm(1);
        for q = 1:k
            mat = mat + lambdaC(j,q)*A{q+1};
            ocena = ocena + abs(lambdaC(j,q))*vecnorm(q+1);
        end
        xr = min_sing_vec(mat,0);
        res = norm(mat*xr);
        % we select eigenvalues that give small residual for the original problem
        if norm(res)<delta*ocena
          lambda = [lambda; lambdaC(j,:)];
          X = [X xr];
      end
    end
else
    error('Unknown method specified in options. Use "compress" or "mep"')
end
    
