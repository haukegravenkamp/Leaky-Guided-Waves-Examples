function [lambda,mu,eta,X1,X2,X3,flag,hist] = threepareigs_si(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,neig,opts)

%THREEPAREIGS_SI  Subspace iteration for a three-parameter eigenvalue problem
% [lambda,mu,X1,X2,X3,flag] = THREEPAREIGS_SI(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,neig,opts)
% returns neig eigenvalues (lambda,mu,eta) with the smallest |eta| of the 
% three-parameter eigenvalue problem. We use suspace compression, TRQI 
% refinement and SVD filtering to expand only with significant new vectors.
%
% A1*x = lambda*B1*x + mu*C1*x + eta*D1*x
% A2*y = lambda*B2*y + mu*C2*y + eta*D2*y
% A3*z = lambda*B3*z + mu*C3*z + eta*D3*z
%
% using subspace iteration with Arnoldi expansion and restart based on 
% selected Ritz vectors on the generalized eigenvalue problem 
% Delta3*w = eta*Delta0*w, 
% where Delta0 and Delta3 are corresponding operator determinants
% Delta0 =   | B1 C1 D1; B2 C2 D2; B3 C3 D3 |
% Delta3 = - | B1 C1 A1; B2 C2 A2; B3 C3 A3 |
%
% Input:
%   - A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3 : matrices, A1, A2, and A3 have to be nonsingular
%   - neig : number of eigenvalues (6)
%   - opts : options (see below)
% 
% Output: 
%   - lambda, mu, eta : eigenvalue parts (eigenvalues are (lambda(j),mu(j),eta(j))
%   - X1, X2, X3 : components of decomposable right eigenvectors 
%     (eigenvector is kron(X1(:,j),kron(X2(:,j),X3(:,j))), such that
%       (A1-lambda(j)*B1-mu(j)*C1-eta(j)*D1)X1(:,j)=0
%       (A2-lambda(j)*B2-mu(j)*C2-eta(j)*D2)X2(:,j)=0
%       (A3-lambda(j)*B3-mu(j)*C3-eta(j)*D3)X3(:,j)=0
%   - flag : convergence (0), no convergence (1)
%   - hist : steps in which eigenvalues were found
%
% Options are (default values in parenthesis):
%   - delta : absolute tolerance (1e-8)*maxinfnorm of the matrices
%   - switcheps : criterion when (norm of the residual) to switch to TRQI refinement for the new eigenpair candidate (1e4*delta)
%   - arnsteps : how many steps we do in the Arnoldi expansion (1)
%   - lowrank : after each iteration we take the leading lowrank Ritz vectors (5)
%   - window : how many Ritz values with minimal |eta| of the projected problem we compute (50)
%   - maxsteps : maximum number of outer steps (500)
%   - maxdetsize : maximum search space (5000), you can enlarge it to e.g. 14000 if you have 16GB RAM
%   - rankeps : in block Krylov expansion we keep only directions from svd with relative sing. values large than rankeps (1e-5) 
%   - etamax : (Inf) we consider only Ritz values with |eta|<=etamax  (default is Inf to consider all)
%   - refine1 : number of TRQI steps in initial Ritz refinement (1)
%   - refine2 : number of TRQI steps to further refine good candidates (3)
%   - showinfo : display Ritz values in each step (2), just eigenvales found (1), nothing (0), default is 1
%   - fp_type: numeric type to use ('single', 'double', or 'mp' (needs MCT),
%     use only if you want to change the default - the superior type of input data
%   - all options for threepareigs
%
% See also: THREEPAREIG, THREEPAREIGS, THREEPAREIGS_JD, TWOPAREIGS_SI,
% DEMO_THREEPAREIGS_SI.

% Reference:
% M.E. Hochstenbach, K. Meerbergen, E. Mengi, B. Plestenjak,
% Subspace methods for 3-parameter eigenvalue problems, arXiv:1802:07386

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 06.09.2015 : no more seed to set random generator
% BP 31.01.2018 : new features - svdfilter, trqi refinement, shrinking
% Modified for MCT

% Last revision: 31.01.2018

narginchk(12,14);
if nargin<13, neig = 6; end
if nargin<14, opts = []; end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3);
end

% Make sure all inputs are of the same numeric type.
if ~isa(A1,class_t), A1 = numeric_t(A1,class_t); end;
if ~isa(B1,class_t), B1 = numeric_t(B1,class_t); end;
if ~isa(C1,class_t), C1 = numeric_t(C1,class_t); end;
if ~isa(D1,class_t), D1 = numeric_t(D1,class_t); end;
if ~isa(A2,class_t), A2 = numeric_t(A2,class_t); end;
if ~isa(B2,class_t), B2 = numeric_t(B2,class_t); end;
if ~isa(C2,class_t), C2 = numeric_t(C2,class_t); end;
if ~isa(D2,class_t), D2 = numeric_t(D2,class_t); end;
if ~isa(A3,class_t), A3 = numeric_t(A3,class_t); end;
if ~isa(B3,class_t), B3 = numeric_t(B3,class_t); end;
if ~isa(C3,class_t), C3 = numeric_t(C3,class_t); end;
if ~isa(D3,class_t), D3 = numeric_t(D3,class_t); end;

maxnorm = max([norm(A1,'inf') norm(B1,'inf') norm(C1,'inf') norm(D1,'inf') ...
   norm(A2,'inf') norm(B2,'inf') norm(C2,'inf') norm(D2,'inf') ...
   norm(A3,'inf') norm(B3,'inf') norm(C3,'inf') norm(D3,'inf')]);

% options for the outer loop
if isfield(opts,'delta'),       delta = opts.delta;              else delta = numeric_t('1e8*eps',class_t)*maxnorm;   end
if isfield(opts,'switcheps'),   switcheps = opts.switcheps;      else switcheps = numeric_t('1e2',class_t)*delta;     end
if isfield(opts,'arnsteps'),    arnsteps = opts.arnsteps;        else arnsteps = 1;           end
if isfield(opts,'maxsteps'),    maxsteps = opts.maxsteps;        else maxsteps = 50;          end
if isfield(opts,'showinfo'),    showinfo = opts.showinfo;        else showinfo = 1;           end
if isfield(opts,'lowrank'),     lowrank = opts.lowrank;          else lowrank = 5;            end
if isfield(opts,'window'),      window = opts.window;            else window = 2*neig;        end
if isfield(opts,'maxdetsize'),  maxdetsize = opts.maxdetsize;    else maxdetsize = 5000;      end
if isfield(opts,'svdfilter'),   svdfilter = opts.svdfilter;      else svdfilter = 1e-5;       end
if isfield(opts,'refine1'),     refine1 = opts.refine1;          else refine1 = 1;            end
if isfield(opts,'refine2'),     refine2 = opts.refine2;          else refine2 = 3;            end
if isfield(opts,'etamax'),      etamax = opts.etamax;            else etamax = Inf;           end
if isfield(opts,'selcrit1'),    selcrit1 = opts.selcrit1;        else selcrit1 = 1e-1;        end
if isfield(opts,'selcrit2'),    selcrit2 = opts.selcrit2;        else selcrit2 = 1e-4;        end
% options for the inner solver threepareigs
if isfield(opts,'tol'),         tol = opts.tol;                  else tol = delta;            end

sparseA = issparse(A1) || issparse(A2) || issparse(A3);

% if all matrices are real we use real search subspaces 
forcereal  = all([isreal(A1) isreal(B1) isreal(C1) isreal(D1) isreal(A2) isreal(B2) isreal(C2) isreal(D2) isreal(A3) isreal(B3) isreal(C3) isreal(D3)]);

lambda = []; mu = []; eta = []; X1 = []; X2 = []; X3 = []; 

optsBlockKrylov.svdfilter = svdfilter;
badTRQI = 0; % we count number of cases when residual does not drop enough after second refine to see if switcheps is too large

% Preparation for Arnoldi expansion and other computations
% If matrices are sparse we solve system in each step, otherwise we precompute the inverse.
if ~sparseA % we explicitly divide by A1 and A2 
   MB1 = A1\B1; MC1 = A1\C1; MD1 = A1\D1;
   MB2 = A2\B2; MC2 = A2\C2; MD2 = A2\D2;
   MB3 = A3\B3; MC3 = A3\C3; MD3 = A3\D3;
end

n1 = size(A1,1); 
n2 = size(A2,1);
n3 = size(A3,1);

% initial random subspaces
U1 = orth(randn(n1,lowrank,class_t));
U2 = orth(randn(n2,lowrank,class_t));
U3 = orth(randn(n3,lowrank,class_t));

% we save left and right eigenvectors of the computed eigenvalues and some
% matrix-vector products from y^'*Delta_0*x products for the selection criteria
EmptyMatrix = numeric_t([],class_t); 
B1XPl = EmptyMatrix; C1XPl = EmptyMatrix; D1XPl = EmptyMatrix; 
B2YPl = EmptyMatrix; C2YPl = EmptyMatrix; D2YPl = EmptyMatrix; 
B3ZPl = EmptyMatrix; C3ZPl = EmptyMatrix; D3ZPl = EmptyMatrix; 
XPl = EmptyMatrix; XPr = EmptyMatrix; 
YPl = EmptyMatrix; YPr = EmptyMatrix; 
ZPl = EmptyMatrix; ZPr = EmptyMatrix; 
InnerProds = EmptyMatrix; 

step = 0;         % number of steps
converged = 0;    % total eigenvalues found so far
oldconverged = 0; % total eigenvalue found before this step 
flag = 0;

% main loop
while (step < maxsteps) && (converged < neig) 
   step = step + 1;
   sizes1 = [size(U1,2) size(U2,2) size(U3,2)];
   % if necessary we expand the search spaces with random vectors
   if sizes1(1)<lowrank,  U1=[U1 randn(n1,lowrank-sizes1(1),class_t)];  end
   if sizes1(2)<lowrank,  U2=[U2 randn(n2,lowrank-sizes1(2),class_t)];  end
   if sizes1(3)<lowrank,  U3=[U3 randn(n3,lowrank-sizes1(3),class_t)];  end
   % Block Arnoldi expansion
   if ~sparseA
        F1 = [MB1*U1 MC1*U1 MD1*U1];
        F2 = [MB2*U2 MC2*U2 MD2*U2];
        F3 = [MB3*U3 MC3*U3 MD3*U3];
        U1n = block_krylov_3p(MB1,MC1,F1,arnsteps,[],optsBlockKrylov);
        U2n = block_krylov_3p(MB2,MC2,F2,arnsteps,[],optsBlockKrylov);
        U3n = block_krylov_3p(MB3,MC3,F3,arnsteps,[],optsBlockKrylov);
   else
        F1 = A1\[B1*U1 C1*U1 D1*U1];
        F2 = A2\[B2*U2 C2*U2 D2*U2];
        F3 = A3\[B3*U3 C3*U3 D3*U3];
        U1n = block_krylov_3p(B1,C1,F1,arnsteps,A1,optsBlockKrylov);
        U2n = block_krylov_3p(B2,C2,F2,arnsteps,A2,optsBlockKrylov);
        U3n = block_krylov_3p(B3,C3,F3,arnsteps,A3,optsBlockKrylov);
   end
   
   sizes2 = [size(U1n,2) size(U2n,2) size(U3n,2)];
   if prod(sizes2)>maxdetsize
       % a proportional reduction of spaces if the projected delta matrices are too large
       vel = shrink(sizes2,maxdetsize);
       U1n = U1n(:,1:vel(1));
       U2n = U2n(:,1:vel(2));
       U3n = U3n(:,1:vel(3));
   end
   sizes3 = [size(U1n,2) size(U2n,2) size(U3n,2)];

   % Projections of initial matrices on the subspace 
   PA1 = A1*U1n;  PB1 = B1*U1n;  PC1 = C1*U1n;  PD1 = D1*U1n;
   PA2 = A2*U2n;  PB2 = B2*U2n;  PC2 = C2*U2n;  PD2 = D2*U2n;
   PA3 = A3*U3n;  PB3 = B3*U3n;  PC3 = C3*U3n;  PD3 = D3*U3n;
   prA1 = U1n'*PA1; prB1 = U1n'*PB1; prC1 = U1n'*PC1; prD1 = U1n'*PD1;
   prA2 = U2n'*PA2; prB2 = U2n'*PB2; prC2 = U2n'*PC2; prD2 = U2n'*PD2;
   prA3 = U3n'*PA3; prB3 = U3n'*PB3; prC3 = U3n'*PC3; prD3 = U3n'*PD3;
   
   % Computation of Ritz values - we compute Ritz values with the smallest % |rho| of the projected 3EP
   [sigma,tau,rho,pXr,pYr,pZr] = threepareigs(prA1,prB1,prC1,prD1,prA2,prB2,prC2,prD2,prA3,prB3,prC3,prD3,window,opts);
   tmpvel = length(sigma);
   res = EmptyMatrix;
   
   % Compute the residuals and refine all candidates using TRQI
   SAXr = U1n*pXr;  SAYr = U2n*pYr;  SAZr = U3n*pZr;
   for j = 1:tmpvel
       if refine1>0
           % we refine Ritz pair using TRQI
           [sigma(j),tau(j),rho(j),SAXr(:,j),SAYr(:,j),SAZr(:,j)] = trqi_3p(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,SAXr(:,j),SAYr(:,j),SAZr(:,j),refine1,tol);
       end
       ABC1 = A1-sigma(j)*B1-tau(j)*C1-rho(j)*D1;
       ABC2 = A2-sigma(j)*B2-tau(j)*C2-rho(j)*D2;
       ABC3 = A3-sigma(j)*B3-tau(j)*C3-rho(j)*D3;
       % norm of the residual
       res(j,1) = norm([ABC1*SAXr(:,j); ABC2*SAYr(:,j); ABC3*SAZr(:,j)], 'fro');	% norm of residual
   end
   minost = min(res);

   % selecton criteria compares Ritz vectors with converged eigenpairs
   if isempty(XPr)
       distance = zeros(tmpvel,1);
       maxprimerjava = zeros(tmpvel,1,class_t);
   else   
       innerprod = (SAXr'*B1XPl).*(SAYr'*C2YPl).*(SAZr'*D3ZPl) + (SAXr'*C1XPl).*(SAYr'*D2YPl).*(SAZr'*B3ZPl) ...
                   + (SAXr'*D1XPl).*(SAYr'*B2YPl).*(SAZr'*C3ZPl) - (SAXr'*B1XPl).*(SAYr'*D2YPl).*(SAZr'*C3ZPl) ...
                   - (SAXr'*C1XPl).*(SAYr'*B2YPl).*(SAZr'*D3ZPl) - (SAXr'*D1XPl).*(SAYr'*C2YPl).*(SAZr'*B3ZPl);
       quotients = abs(innerprod)./(ones(tmpvel,1,class_t)*InnerProds);
       maxprimerjava = max(quotients,[],2);
   end
   newstepindex = [];
   
   flag = zeros(tmpvel,1); % status of Ritz pairs
   
   for j = 1:tmpvel
       if (maxprimerjava(j)<selcrit1) % selection criteria is the first check
            flag(j,1) = 2; % Ritz pair satisfied tzhe selection criteria
            if res(j)<switcheps && ( abs(rho(j))<etamax ) % this is a good candidate
                flag(j,1) = 3; % this is a good candidate for an eigenpair
                % further refine good candidates using additional TRQI steps
                if refine2>0
                    [la,ua,ea,Xrapr,Yrapr,Zrapr] = trqi_3p(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,SAXr(:,j),SAYr(:,j),SAZr(:,j),refine2,tol);
                    SAXr(:,j) = Xrapr; 
                    SAYr(:,j) = Yrapr;
                    SAZr(:,j) = Zrapr;
                else
                    la = sigma(j); ua = tau(j); ea = rho(j);
                    Xrapr = SAXr(:,j); 
                    Yrapr = SAYr(:,j);
                    Zrapr = SAZr(:,j);
                end
                ABC1 = A1-la*B1-ua*C1-ea*D1;
                ABC2 = A2-la*B2-ua*C2-ea*D2;
                ABC3 = A3-la*B3-ua*C3-ea*D3;
                resnorm = norm([ABC1*Xrapr; ABC2*Yrapr; ABC3*Zrapr], 'fro'); 	% norm of the residual
                
                if (refine2>0) && (resnorm>=delta)
                    badTRQI = badTRQI + 1;
                end

                if isempty(XPr)
                    quotients = 0;
                else   
                    innerprod = (Xrapr'*B1XPl).*(Yrapr'*C2YPl).*(Zrapr'*D3ZPl) + (Xrapr'*C1XPl).*(Yrapr'*D2YPl).*(Zrapr'*B3ZPl) ...
                              + (Xrapr'*D1XPl).*(Yrapr'*B2YPl).*(Zrapr'*C3ZPl) - (Xrapr'*B1XPl).*(Yrapr'*D2YPl).*(Zrapr'*C3ZPl) ...
                              - (Xrapr'*C1XPl).*(Yrapr'*B2YPl).*(Zrapr'*D3ZPl) - (Xrapr'*D1XPl).*(Yrapr'*C2YPl).*(Zrapr'*B3ZPl);
                    quotients = max(abs(innerprod)./InnerProds);
                end
            
                % we check the residual and the selection criteria again (as both might change because of TRQI)
                if (resnorm<delta) && (quotients<selcrit2) 
                    flag(j,1) = 4; % converged eigenvalue
                    XPr = [XPr Xrapr];
                    YPr = [YPr Yrapr];
                    ZPr = [ZPr Zrapr];
                    X1 = [X1 Xrapr];
                    X2 = [X2 Yrapr];
                    X3 = [X3 Zrapr];
                    % left eigenvectors are computed for the selection criteria
                    [tilda, Xlapr] = min_sing_vec(ABC1,1);
                    [tilda, Ylapr] = min_sing_vec(ABC2,1);
                    [tilda, Zlapr] = min_sing_vec(ABC3,1);
                    B1Xlapr = B1*Xlapr; C1Xlapr = C1*Xlapr; D1Xlapr = D1*Xlapr;
                    B2Ylapr = B2*Ylapr; C2Ylapr = C2*Ylapr; D2Ylapr = D2*Ylapr;
                    B3Zlapr = B3*Zlapr; C3Zlapr = C3*Zlapr; D3Zlapr = D3*Zlapr;
                    XPl = [XPl Xlapr];
                    YPl = [YPl Ylapr];
                    ZPl = [ZPl Zlapr];
                    B1XPl = B1'*XPl; C1XPl = C1'*XPl; D1XPl = D1'*XPl;
                    B2YPl = B2'*YPl; C2YPl = C2'*YPl; D2YPl = D2'*YPl;
                    B3ZPl = B3'*ZPl; C3ZPl = C3'*ZPl; D3ZPl = D3'*ZPl;
                    innerpr = (Xrapr'*B1Xlapr).*(Yrapr'*C2Ylapr).*(Zrapr'*D3Zlapr) + (Xrapr'*C1Xlapr).*(Yrapr'*D2Ylapr).*(Zrapr'*B3Zlapr) ...
                        + (Xrapr'*D1Xlapr).*(Yrapr'*B2Ylapr).*(Zrapr'*C3Zlapr) - (Xrapr'*B1Xlapr).*(Yrapr'*D2Ylapr).*(Zrapr'*C3Zlapr) ...
                        - (Xrapr'*C1Xlapr).*(Yrapr'*B2Ylapr).*(Zrapr'*D3Zlapr) - (Xrapr'*D1Xlapr).*(Yrapr'*C2Ylapr).*(Zrapr'*B3Zlapr);
                    InnerProds = [InnerProds abs(innerpr)];
                    converged = converged + 1;
                    lambda = [lambda; la];  
                    mu = [mu; ua];  
                    eta = [eta; ea];               
                    disp(sprintf('Eig (%2d): lambda: %11.4e%+11.4ei, mu: %11.4e%+11.4ei, eta: %11.4e%+11.4ei, step %4d, inres: %5.1e, res: %5.1e, selcrit: %5.1e',converged,real(la),imag(la),real(ua),imag(ua),real(ea),imag(ea),step,res(j),resnorm,quotients))
                else
                    if (quotients<selcrit2) 
                        newstepindex = [newstepindex j]; % residual was not small enough for convergence, selection criteria is still satisfied
                    else
                        flag(j,1) = 5; % method converged to an already computed eigenvalue (repeated)
                        res(j) = Inf;
                    end
                end

            else
               % we take this for new direction    
               newstepindex = [newstepindex j];
            end
            
       else
          flag(j,1) = 6; % rejected by the selection criteria
          res(j) = Inf;
       end
       if (converged == neig) 
           break
       end
   end
    if length(newstepindex)>lowrank
        newstepindex=newstepindex(1:lowrank);
    end
        
    if showinfo == 2
        sweep = [sigma(1:tmpvel) tau(1:tmpvel) rho(1:tmpvel) res flag] 
    end
    rejected = sum(isinf(res(:,1)));
    repeated = sum(flag==5);
    if showinfo
       fprintf('Step %3d, from: %d x %d x %d -> b. Krylov: %2d x %2d x %2d (%5d) -> red.: %2d x %2d x %2d (%4d), found %2d, new %2d, rejected %2d, repeated %2d, minost: %4.1e\n',...
          step,sizes1(1),sizes1(2),sizes1(3),sizes2(1),sizes2(2),sizes2(3),prod(sizes2),sizes3(1),sizes3(2),sizes3(3),prod(sizes3),converged,converged-oldconverged,...
          rejected,repeated,minost)
    end
    if converged < neig
       % we form new subspace out of current approximations
       V1next = SAXr(:,newstepindex); 
       V2next = SAYr(:,newstepindex); 
       V3next = SAZr(:,newstepindex); 
       if forcereal
           % for real matrices we use real subspaces
           % this enables us to use larger subspaces as real problems require less memory
           U1 = orth([real(V1next) imag(V1next)]);
           U2 = orth([real(V2next) imag(V2next)]);
           U3 = orth([real(V3next) imag(V3next)]);
       else
           U1 = orth(V1next);
           U2 = orth(V2next);
           U3 = orth(V3next);
       end
       oldconverged = converged;
    end    
end

if converged < neig
    fprintf('Method has not converged, returning %d eigenvalues from last step\n',converged);
    flag = 1;
end

hist.steps = step;
hist.badTRQI = badTRQI;

end %% threepareigs_si

% ------------------------------------------------------------------------
% Auxiliary function shrinks the subspaces by deleting a proportional 
% number of last columns in each basis
% ------------------------------------------------------------------------

function b = shrink(a,n)

[v,p] = sort(a);
vel = prod(v);
if vel <= n
    b = a;
else
    fac = (vel/n)^(1/3);
    w(1) = round(v(1)/fac);
    fac = (v(2)*v(3)*w(1)/n)^(1/2);
    w(2) = round(v(2)/fac);
    w(3) = floor(n/(w(1)*w(2)));
    [tmp,invp] = sort(p);
    b = w(invp);
end
end
