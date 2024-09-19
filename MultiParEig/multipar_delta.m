function Delta = multipar_delta(A,indices)

%MULTIPAR_DELTA   Delta matrices for a multiparameter eigenvalue problem
%
% Delta = MULTIPAR_DELTA(A) returns set of k+1 Delta matrices, which are 
% operator determinants related to the multiparameter eigenvalue problem
%
% A{1,1} x1 = lambda(1) A{1,2} x1 +  ... + lambda(k) A{1,k+1} x1 
% A{2,1} x2 = lambda(1) A{2,2} x2 +  ... + lambda(k) A{2,k+1} x2 
% ...
% A{k,1} xk = lambda(1) A{k,2} xk +  ... + lambda(k) A{k,k+1} xk 
%
% See also: TWOPAR_DELTA, THREEPAR_DELTA
%
% If you want to compute just some of the matrices Delta, use 
% Delta = MULTIPAR_DELTA(A,indices), where indices is a vector with indices
% 0 to k of Delta matrices that you want to compute.

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% First revision: 3.7.2022
% BP 08.11.2023 : added direct computations for k=2,3,4

    k = length(A) - 1; % number of parameters 
    if nargin<2
        indices = 0:k;
    end

    Delta = cell(1,k+1);
    if k==2
        [D1,D2,D3] = twopar_delta(A{1,1},A{1,2},A{1,3},A{2,1},A{2,2},A{2,3},indices);
        Delta = cell(1,3);
        Delta{1} = D1;
        Delta{2} = D2;
        Delta{3} = D3;
    elseif k==3
        [D1,D2,D3,D4] = threepar_delta(A{1,1},A{1,2},A{1,3},A{1,4},A{2,1},A{2,2},A{2,3},A{2,4},A{3,1},A{3,2},A{3,3},A{3,4},indices);
        Delta = cell(1,4);
        Delta{1} = D1;
        Delta{2} = D2;
        Delta{3} = D3;
        Delta{4} = D4;
    elseif k==4
        [D1,D2,D3,D4,D5] = fourpar_delta(A{1,1},A{1,2},A{1,3},A{1,4},A{1,5},A{2,1},A{2,2},A{2,3},A{2,4},A{2,5},A{3,1},A{3,2},A{3,3},A{3,4},A{3,5},A{4,1},A{4,2},A{4,3},A{4,4},A{4,5},indices);
        Delta = cell(1,5);
        Delta{1} = D1;
        Delta{2} = D2;
        Delta{3} = D3;
        Delta{4} = D4;
        Delta{5} = D5;
    else
        % computation of Delta matrices
        for j = 1:length(indices)
            r = indices(j);
            Delta{r+1} = krondet(A, r);
        end
    end

end

%------------------------------------------------------------------------

function krondelta = krondet(A,index)
% recursive computation of Delta matrices - operator determinants

    if nargin == 2
        % we replace index-th column with the first (if index>0)
        n = length(A);
        if index==0
            ind = 2:n;
        else
            ind = [2:index 1 (index+2):n];
        end
        A = A(:, ind);
    end
    n = length(A);
    if n == 1
        krondelta = A{1,1}; 
        return;
    end
    deter = []; 
    sgn = 1;
    for k = 1:n
        indj = [1:(k-1), (k+1):n];
        indi = 2:n;   
        if isempty(deter)
            deter = sgn*kron( A{1, k}, krondet(A(indi, indj)) );
        else
            deter = deter + sgn*kron( A{1, k}, krondet(A(indi, indj)) );
        end
        sgn = -sgn;
    end
    
    krondelta = deter;
    
end

%------------------------------------------------------------------------
function [Delta0,Delta1,Delta2] = twopar_delta(A1,B1,C1,A2,B2,C2,indices)
% direct formulas for 2 x 2 determinants

    if any(indices==0)
        Delta0 = kron(B1,C2) - kron(C1,B2);
    else
        Delta0 = [];
    end
    
    if any(indices==1)
        Delta1 = kron(A1,C2) - kron(C1,A2);
    else
        Delta1 = [];
    end
    
    if any(indices==2)
        Delta2 = kron(B1,A2) - kron(A1,B2);
    else
        Delta2 = [];
    end

end
%------------------------------------------------------------------------
function [Delta0,Delta1,Delta2,Delta3] = threepar_delta(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,indices)
% direct formulas for 3 x 3 determinants

    if any(indices==0)
        Delta0 = kron(B1,kron(C2,D3)-kron(D2,C3)) - kron(C1,kron(B2,D3)-kron(D2,B3)) + kron(D1,kron(B2,C3) - kron(C2,B3));
    else
        Delta0 = [];
    end
    
    if any(indices==1)
        Delta1 = kron(A1,kron(C2,D3)-kron(D2,C3)) - kron(C1,kron(A2,D3)-kron(D2,A3)) + kron(D1,kron(A2,C3) - kron(C2,A3));
    else
        Delta1 = [];
    end
    
    if any(indices==2)
        Delta2 = kron(B1,kron(A2,D3)-kron(D2,A3)) - kron(A1,kron(B2,D3)-kron(D2,B3)) + kron(D1,kron(B2,A3) - kron(A2,B3));
    else
        Delta2 = [];
    end
    
    if any(indices==3)
        Delta3 = kron(B1,kron(C2,A3)-kron(A2,C3)) - kron(C1,kron(B2,A3)-kron(A2,B3)) + kron(A1,kron(B2,C3) - kron(C2,B3));
    else
        Delta3 = [];
    end

end
%------------------------------------------------------------------------
function [Delta0,Delta1,Delta2,Delta3,Delta4] = fourpar_delta(A1,B1,C1,D1,E1,A2,B2,C2,D2,E2,A3,B3,C3,D3,E3,A4,B4,C4,D4,E4,indices)
% direct formulas for 4 x 4 determinants

    if any(indices==0)
        Delta0 = kron(B1,kron(C2,kron(D3,E4)-kron(E3,D4)) - kron(D2,kron(C3,E4)-kron(E3,C4)) + kron(E2,kron(C3,D4) - kron(D3,C4))) - kron(C1,kron(B2,kron(D3,E4)-kron(E3,D4)) - kron(D2,kron(B3,E4)-kron(E3,B4)) + kron(E2,kron(B3,D4) - kron(D3,B4))) + kron(D1,kron(B2,kron(C3,E4)-kron(E3,C4)) - kron(C2,kron(B3,E4)-kron(E3,B4)) + kron(E2,kron(B3,C4) - kron(C3,B4))) - kron(E1,kron(B2,kron(C3,D4)-kron(D3,C4)) - kron(C2,kron(B3,D4)-kron(D3,B4)) + kron(D2,kron(B3,C4) - kron(C3,B4))); 
    else
        Delta0 = [];
    end
    
    if any(indices==1)
        Delta1 = kron(A1,kron(C2,kron(D3,E4)-kron(E3,D4)) - kron(D2,kron(C3,E4)-kron(E3,C4)) + kron(E2,kron(C3,D4) - kron(D3,C4))) - kron(C1,kron(A2,kron(D3,E4)-kron(E3,D4)) - kron(D2,kron(A3,E4)-kron(E3,A4)) + kron(E2,kron(A3,D4) - kron(D3,A4))) + kron(D1,kron(A2,kron(C3,E4)-kron(E3,C4)) - kron(C2,kron(A3,E4)-kron(E3,A4)) + kron(E2,kron(A3,C4) - kron(C3,A4))) - kron(E1,kron(A2,kron(C3,D4)-kron(D3,C4)) - kron(C2,kron(A3,D4)-kron(D3,A4)) + kron(D2,kron(A3,C4) - kron(C3,A4))); 
    else
        Delta1 = [];
    end
    
    if any(indices==2)
        Delta2 = kron(B1,kron(A2,kron(D3,E4)-kron(E3,D4)) - kron(D2,kron(A3,E4)-kron(E3,A4)) + kron(E2,kron(A3,D4) - kron(D3,A4))) - kron(A1,kron(B2,kron(D3,E4)-kron(E3,D4)) - kron(D2,kron(B3,E4)-kron(E3,B4)) + kron(E2,kron(B3,D4) - kron(D3,B4))) + kron(D1,kron(B2,kron(A3,E4)-kron(E3,A4)) - kron(A2,kron(B3,E4)-kron(E3,B4)) + kron(E2,kron(B3,A4) - kron(A3,B4))) - kron(E1,kron(B2,kron(A3,D4)-kron(D3,A4)) - kron(A2,kron(B3,D4)-kron(D3,B4)) + kron(D2,kron(B3,A4) - kron(A3,B4))); 
    else
        Delta2 = [];
    end
    
    if any(indices==3)
        Delta3 = kron(B1,kron(C2,kron(A3,E4)-kron(E3,A4)) - kron(A2,kron(C3,E4)-kron(E3,C4)) + kron(E2,kron(C3,A4) - kron(A3,C4))) - kron(C1,kron(B2,kron(A3,E4)-kron(E3,A4)) - kron(A2,kron(B3,E4)-kron(E3,B4)) + kron(E2,kron(B3,A4) - kron(A3,B4))) + kron(A1,kron(B2,kron(C3,E4)-kron(E3,C4)) - kron(C2,kron(B3,E4)-kron(E3,B4)) + kron(E2,kron(B3,C4) - kron(C3,B4))) - kron(E1,kron(B2,kron(C3,A4)-kron(A3,C4)) - kron(C2,kron(B3,A4)-kron(A3,B4)) + kron(A2,kron(B3,C4) - kron(C3,B4))); 
    else
        Delta3 = [];
    end
    
    if any(indices==4)
        Delta4 = kron(B1,kron(C2,kron(D3,A4)-kron(A3,D4)) - kron(D2,kron(C3,A4)-kron(A3,C4)) + kron(A2,kron(C3,D4) - kron(D3,C4))) - kron(C1,kron(B2,kron(D3,A4)-kron(A3,D4)) - kron(D2,kron(B3,A4)-kron(A3,B4)) + kron(A2,kron(B3,D4) - kron(D3,B4))) + kron(D1,kron(B2,kron(C3,A4)-kron(A3,C4)) - kron(C2,kron(B3,A4)-kron(A3,B4)) + kron(A2,kron(B3,C4) - kron(C3,B4))) - kron(A1,kron(B2,kron(C3,D4)-kron(D3,C4)) - kron(C2,kron(B3,D4)-kron(D3,B4)) + kron(D2,kron(B3,C4) - kron(C3,B4))); 
    else
        Delta4 = [];
    end

end