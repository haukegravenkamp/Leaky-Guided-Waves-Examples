function [lambda,mu,eta,XR,YR,ZR,XL,YL,ZL,res,hist] = threepareigs_jd(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,neig,opts)

%THREEPAREIGS_JD   Jacobi-Davidson method for a three-parameter eigenvalue problem
%
% [lambda,mu,eta,XR,YR,ZR,XL,YL,ZL,res] = THREEPAREIGS_JD(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,neig,opts)
% returns neig eigenvalues of the three-parameter eigenvalue problem
%
% A1*x = lambda*B1*x + mu*C1*x + eta*D1*x
% A2*y = lambda*B2*y + mu*C2*y + eta*D2*y
% A3*z = lambda*B3*z + mu*C3*z + eta*D3*z
%
% using the Jacobi-Davidson method.
%
% Input:
%   - A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3 : matrices
%   - neig : number of eigenvalues (1)             
%   - opts : options (see below)
%
% Output:
%   - lambda, mu, eta : eigenvalue parts (eigenvalues are (lambda(j),mu(j),eta(j))
%   - XR, YR, ZR : components of decomposable right eigenvectors 
%     (eigenvector is kron(XR(:,j),kron(YR(:,j),ZR(:,j))), such that
%       (A1-lambda(j)*B1-mu(j)*C1-eta(j)*D1)*XR(:,j)=0
%       (A2-lambda(j)*B2-mu(j)*C2-eta(j)*D2)*YR(:,j)=0
%       (A3-lambda(j)*B3-mu(j)*C3-eta(j)*D3)*ZR(:,j)=0
%   - XL, YL, ZL : components of decomposable left eigenvectors 
%     (eigenvector is kron(XL(:,j),kron(YL(:,j),ZL(:,j))), such that
%       (A1-lambda(j)*B1-mu(j)*C1-eta(j)*D1)'*XL(:,j)=0
%       (A2-lambda(j)*B2-mu(j)*C2-eta(j)*D2)'*YL(:,j)=0
%       (A3-lambda(j)*B3-mu(j)*C3-eta(j)*D3)'*ZL(:,j)=0
%   - res : norms of residuals for all J-D steps, try plot(log10(res))
%   - hist : additional information for developers
%
% Options are (default values in parenthesis):
%   - delta : convergence criterion (norm of the residual) for the outer iteration (eps*max inf norm of the matrices)
%   - switcheps : criterion when (norm of the residual) to switch to TRQI refinement for the new eigenpair candidate (1e4*delta)
%   - X0r,Y0r,Z0r : initial search space (rand(n1,1), rand(n2,1), rand(n3,1))
%   - minsize : dimension of search spaces after restart (5)
%   - maxsize : maximum dimension of search spaces before restart (10)
%   - maxsteps : maximum number of outer iteration steps (100)
%   - innersteps : number of GMRES steps for the correction equation (0, set to -1 to solve the correction equation exactly)
%   - innertol : tolerance in GMRES for the correction equation (0)
%   - window : number of Ritz values with the smallest |mu| that we compute, set to 0 (default) for all Ritz values
%   - extraction : extraction method, choices are: 
%        'maxdist' : the maximum distance from the target,
%        'mindist' : the minimum distance from the target (default),
%        'minres'  : the smallest residual,
%        'maxlambda': the eigenvalue with lambda with the largest real part
%        'mineta' : the eigenvalue with the smallest |eta|
%   - target : target for the eigenvalues ([0 0 0])
%   - reschange : switch to minimum residual extraction when residual norm is small - (10^(-2.5))
%   - XPr,XPl,YPl,YPr,ZPr,ZPl  - prohibited directions (previously computed left and right eigenvectors) - ([],[],[],[],[],[])
%   - M1, M2, M3 : left preconditioners -  can also be a function_handle that returns M1*x1, M2*x2, or M3*x3
%   - harmonic : set to 1 to use harmonic instead of Ritz values (0) - use this for interior eigenvalue
%   - forcereal : set to 1 if you know that eigenvectors and eigenvalues are real, default is 0
%   - refine : (3) number of TRQI steps to refine an eigenpair candidate
%   - solveropts : options for threepareigs of threepareigs ([])
%   - showinfo : display nothing (default), 1: all values, 2: just eigenvalues
%   - fp_type: numeric type to use ('single', 'double', or 'mp' (needs MCT),
%     use only if you want to change the default - the superior type of input data
%
% See also: THREEPAREIG, THREEPAREIGS, THREEPAREIGS_SI, TWOPAREIGS_JD

% References:
%
%   1) M.E. Hochstenbach, B. Plestenjak,
%      A Jacobi-Davidson type method for a right definite two-parameter eigenvalue problem,
%      SIAM J. Matrix Anal. Appl. 24 (2002), 392-410.
%
%   2) M.E. Hochstenbach, B. Plestenjak, T. Kosir,
%      A Jacobi-Davidson type method for the two-parameter eigenvalue problem,
%      SIAM J. Matrix Anal. Appl. 26 (2005), 477-497.
%
%   3) M.E. Hochstenbach, B. Plestenjak,
%      Harmonic Rayleigh-Ritz extraction for the multiparameter eigenvalue problem,
%      Electron. Trans. Numer. Anal. 29 (2008) 81-96.
%
%   4) M.E. Hochstenbach, K. Meerbergen, E. Mengi, B. Plestenjak,
%      Subspace methods for 3-parameter eigenvalue problems,
%      arXiv:1802:07386

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 06.11.2016 : Switch to TRQI refinement and some small changes
% BP Fixed bug in the preconditioned first order correction 
% BP 31.01.208 : new features - exact solution, trqi refinement, new restart
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

n1 = size(A1,1); 
n2 = size(A2,1);
n3 = size(A3,1);

EmptyMatrix = numeric_t([],class_t);

if isfield(opts,'delta')
    delta = getfield(opts,'delta');            
else
    maxnor = max([norm(A1,'inf') norm(B1,'inf') norm(C1,'inf') norm(D1,'inf') ...
                norm(A2,'inf') norm(B2,'inf') norm(C2,'inf') norm(D2,'inf') ...
                norm(A3,'inf') norm(B3,'inf') norm(C3,'inf') norm(D3,'inf')]);
    delta = numeric_t('10*eps',class_t)*maxnor;            
end
if isfield(opts,'switcheps'),   switcheps = getfield(opts,'switcheps');    else switcheps = numeric_t('1e4',class_t)*delta;   end
if isfield(opts,'minsize'),     minsize = getfield(opts,'minsize');        else minsize = 5;             end
if isfield(opts,'maxsize'),     maxsize = getfield(opts,'maxsize');        else maxsize = 10;            end
if isfield(opts,'maxsteps'),    maxsteps = getfield(opts,'maxsteps');      else maxsteps = 100;          end
if isfield(opts,'showinfo'),    showinfo = getfield(opts,'showinfo');      else showinfo = 1;            end
if isfield(opts,'target'),      target = getfield(opts,'target');          else target = [0 0 0];        end
if isfield(opts,'X0r'),         X0r = getfield(opts,'X0r');                else X0r = rand(n1,1,class_t); end
if isfield(opts,'Y0r'),         Y0r = getfield(opts,'Y0r');                else Y0r = rand(n2,1,class_t); end
if isfield(opts,'Z0r'),         Z0r = getfield(opts,'Z0r');                else Z0r = rand(n3,1,class_t); end
if isfield(opts,'innersteps'),  innersteps = getfield(opts,'innersteps');  else innersteps = 0;          end
if isfield(opts,'innertol'),    innertol = getfield(opts,'innertol');      else innertol = 0;            end
if isfield(opts,'window'),      window = getfield(opts,'window');          else window = 0;              end
if isfield(opts,'extraction'),  extraction = getfield(opts,'extraction');  else extraction = 'mindist';  end
if isfield(opts,'reschange'),   reschange = getfield(opts,'reschange');    else reschange = 10^(-2.5);   end
if isfield(opts,'XPr'),         XPr = getfield(opts,'XPr');                else XPr = EmptyMatrix;       end
if isfield(opts,'XPl'),         XPl = getfield(opts,'XPl');                else XPl = EmptyMatrix;       end
if isfield(opts,'YPr'),         YPr = getfield(opts,'YPr');                else YPr = EmptyMatrix;       end
if isfield(opts,'YPl'),         YPl = getfield(opts,'YPl');                else YPl = EmptyMatrix;       end
if isfield(opts,'ZPr'),         ZPr = getfield(opts,'ZPr');                else ZPr = EmptyMatrix;       end
if isfield(opts,'ZPl'),         ZPl = getfield(opts,'ZPl');                else ZPl = EmptyMatrix;       end
if isfield(opts,'M1'),          M1 = getfield(opts,'M1');                  else M1 = [];                 end
if isfield(opts,'M2'),          M2 = getfield(opts,'M2');                  else M2 = [];                 end
if isfield(opts,'M3'),          M3 = getfield(opts,'M3');                  else M3 = [];                 end
if isfield(opts,'harmonic'),    harmonic = getfield(opts,'harmonic');      else harmonic = 0;            end
if isfield(opts,'forcereal'),   forcereal = getfield(opts,'forcereal');    else forcereal = 0;           end
if isfield(opts,'refine'),      refine = getfield(opts,'refine');          else refine = 4;              end
if isfield(opts,'solveropts'),  solveropts = getfield(opts,'solveropts');  else solveropts = [];         end
if isfield(opts,'selcrit1'),    selcrit1 = getfield(opts,'selcrit1');      else selcrit1 = 1e-1;         end
if isfield(opts,'selcrit2'),    selcrit2 = getfield(opts,'selcrit2');      else selcrit2 = 1e-4;         end

if ~isfield(solveropts,'usesparse'),  solveropts.usesparse = 0; end;

EmptyMatrix = numeric_t([],class_t);

if ~isempty(XPl)
   B1XPl = B1'*XPl; C1XPl = C1'*XPl; D1XPl = D1'*XPl;
   B2YPl = B2'*YPl; C2YPl = C2'*YPl; D2YPl = D2'*YPl;
   B3YPl = B3'*ZPl; C3ZPl = C3'*ZPl; D3ZPl = D3'*ZPl;
else
   B1XPl = EmptyMatrix; C1XPl = EmptyMatrix; D1XPl = EmptyMatrix;
   B2YPl = EmptyMatrix; C2YPl = EmptyMatrix; D2YPl = EmptyMatrix;
   B3ZPl = EmptyMatrix; C3ZPl = EmptyMatrix; D3ZPl = EmptyMatrix;
end
% Initial search spaces and other initialization
U1 = rgs(X0r); 
U2 = rgs(Y0r); 
U3 = rgs(Z0r); 
RestX = EmptyMatrix; RestY = EmptyMatrix; RestZ = EmptyMatrix;

initialextraction = extraction;
lambda = EmptyMatrix; mu = EmptyMatrix; eta = EmptyMatrix;
conv = 0; % no convergence yet
step = 1; % number of steps
maxsteps = maxsteps + 1; % as we start counting with 1 

lastconv = 0; % step with the last convergence

AU1 = A1*U1;  BU1 = B1*U1;  CU1 = C1*U1;  DU1 = D1*U1;
AU2 = A2*U2;  BU2 = B2*U2;  CU2 = C2*U2;  DU2 = D2*U2;
AU3 = A3*U3;  BU3 = B3*U3;  CU3 = C3*U3;  DU3 = D3*U3;

if harmonic
    HABC1 = A1-target(1)*B1-target(2)*C1-target(3)*D1;
    HABC2 = A2-target(1)*B2-target(2)*C2-target(3)*D2;
    HABC3 = A3-target(1)*B3-target(2)*C3-target(3)*D3;
    HABCU1 = HABC1*U1;
    HABCU2 = HABC2*U2;
    HABCU3 = HABC3*U3;
    Q1 = rgs(HABCU1); 
    Q2 = rgs(HABCU2); 
    Q3 = rgs(HABCU3); 
end    

% selection criteria for Delta0 orthogonality
InnerProds = EmptyMatrix;
if ~isempty(XPr) 
   for j = 1:size(Xpr,2)
        innerpr = (XPr(:,j)'*B1XPl(:,j)).*(YPr(:,j)'*C2YPl(:,j)).*(ZPr(:,j)'*D3ZPl(:,j)) + (XPr(:,j)'*C1XPl(:,j)).*(YPr(:,j)'*D2YPl(:,j)).*(ZPr(:,j)'*B3ZPl(:,j)) ...
             + (XPr(:,j)'*D1XPl(:,j)).*(YPr(:,j)'*B2YPl(:,j)).*(ZPr(:,j)'*C3ZPl(:,j)) - (XPr(:,j)'*B1XPl(:,j)).*(YPr(:,j)'*D2YPl(:,j)).*(ZPr(:,j)'*C3ZPl(:,j)) ...
             - (XPr(:,j)'*C1XPl(:,j)).*(YPr(:,j)'*B2YPl(:,j)).*(ZPr(:,j)'*D3ZPl(:,j)) - (XPr(:,j)'*D1XPl(:,j)).*(YPr(:,j)'*C2YPl(:,j)).*(ZPr(:,j)'*B3ZPl(:,j));
        InnerProds = [InnerProds abs(innerpr)];
   end
   minInnerProd = min(InnerProds);
else
   minInnerProd = Inf; 
end

res = EmptyMatrix;

hist.missedTRQI = 0;
hist.repetead = 0;

% Start of the method (loop of restart cycles)
while (conv<neig) && (step<maxsteps) 
      
  while (size(U1,2)<=maxsize) && (conv<neig) && (step<maxsteps) 
      
    % Step 1: We find appropriate (harmonic) Ritz value in search space
    % --------------------------------------------------------------------------------------
    % Solution of smaller problem in the search base
    
    if harmonic
        SA1 = Q1'*HABCU1;  SB1 = Q1'*BU1;  SC1 = Q1'*CU1;  SD1 = Q1'*DU1;
        SA2 = Q2'*HABCU2;  SB2 = Q2'*BU2;  SC2 = Q2'*CU2;  SD2 = Q2'*DU2;
        SA3 = Q3'*HABCU3;  SB3 = Q3'*BU3;  SC3 = Q3'*CU3;  SD3 = Q3'*DU3;
    else
        SA1 = U1'*AU1;  SB1 = U1'*BU1;  SC1 = U1'*CU1;  SD1 = U1'*DU1;
        SA2 = U2'*AU2;  SB2 = U2'*BU2;  SC2 = U2'*CU2;  SD2 = U2'*DU2;
        SA3 = U3'*AU3;  SB3 = U3'*BU3;  SC3 = U3'*CU3;  SD3 = U3'*DU3;
    end
    
    if window>0
        [Al,Au,Ae,AXr,AYr,AZr] = threepareigs(SA1,SB1,SC1,SD1,SA2,SB2,SC2,SD2,SA3,SB3,SC3,SD3,window,solveropts);
    else
        [Al,Au,Ae,AXr,AYr,AZr] = threepareig(SA1,SB1,SC1,SD1,SA2,SB2,SC2,SD2,SA3,SB3,SC3,SD3,solveropts);
    end
    
    tmpvel = length(Al);

    if harmonic
        Al = Al + target(1);
        Au = Au + target(2);
        Ae = Ae + target(3);
    end

    noexpansion = 1;   
    firstcand = 1;
    
    % in a loop we can extract more than one eigenpair from the subspace
    while (conv<neig) && noexpansion 
        
       switch extraction
           case 'maxdist'  % Ritz value farthest to the target
               ritznorm = abs(Al-target(1)).^2 + abs(Au-target(2)).^2 + abs(Ae-target(3)).^2;
               [tilda,order] = sort(-ritznorm);
           case 'mineta'  % Ritz value with the smallest |eta|
               [tilda,order] = sort(abs(Ae));
           case 'mindist'  % Ritz value closest to the target 
               ritznorm = abs(Al-target(1)).^2 + abs(Au-target(2)).^2 + abs(Ae-target(3)).^2;
               [tilda,order] = sort(ritznorm);
           case 'maxlambda'  % Ritz value with the largest real part
               [tilda,order] = sort(-Al);
           case 'minres'  % Ritz pair with the smallest residual               
               rn = zeros(length(Al),1);
               for kk = 1:length(Al)
                   ABCU1 = AU1-Al(kk)*BU1-Au(kk)*CU1-Ae(kk)*DU1;
                   ABCU2 = AU2-Al(kk)*BU2-Au(kk)*CU2-Ae(kk)*DU2;
                   ABCU3 = AU3-Al(kk)*BU3-Au(kk)*CU3-Ae(kk)*DU3;
                   r1r = ABCU1*AXr(:,kk); % residual r1 right
                   r2r = ABCU2*AYr(:,kk); % residual r2 right
                   r3r = ABCU3*AZr(:,kk); % residual r3 right
                   rn(kk,1) = norm([r1r;r2r;r3r], 'fro'); % norm of residual
               end   
               [tilda,order] = sort(rn);
           otherwise
               error('Unknown extraction option')
        end
         
        % we select the Ritz vector that is Delta0 orthogonal to computed eigenpairs
        if isempty(XPr)
            selcriteria = zeros(length(Al),1);
        else   
            SAXr = U1*AXr;  SAYr = U2*AYr;  SAZr = U3*AZr;
            innerprod =   (SAXr'*B1XPl).*(SAYr'*C2YPl).*(SAZr'*D3ZPl) + (SAXr'*C1XPl).*(SAYr'*D2YPl).*(SAZr'*B3ZPl) ...
                        + (SAXr'*D1XPl).*(SAYr'*B2YPl).*(SAZr'*C3ZPl) - (SAXr'*B1XPl).*(SAYr'*D2YPl).*(SAZr'*C3ZPl) ...
                        - (SAXr'*C1XPl).*(SAYr'*B2YPl).*(SAZr'*D3ZPl) - (SAXr'*D1XPl).*(SAYr'*C2YPl).*(SAZr'*B3ZPl);
            primerjava = abs(innerprod)./(ones(tmpvel,1)*InnerProds);
            selcriteria = max(primerjava,[],2);
        end
        
        izb = firstcand;
        while selcriteria(order(izb))>selcrit1
            izb = izb+1;
            if izb>length(order);  break;  end
        end
        skipped = izb-1;
   
        if izb>length(order)  
            % all directions are not enough "Delta0 orthogonal", so we pick the one with the minimal distance    
            [tilda,pos] = min(selcriteria);
            izb = 1;
            expandrandom = 1;  % safety flag that prevents taking already computed eigenvalue
        else
            pos = order(izb);
            expandrandom = 0;
        end;
        
        dddd = selcriteria(pos); 
        la = Al(pos);  ua = Au(pos);  ea = Ae(pos);
        Xrapr = U1*AXr(:,pos);  Yrapr = U2*AYr(:,pos); Zrapr = U3*AZr(:,pos);  % best Ritz vectors
        
        if harmonic
            % TRQ improvement
            A1Xr = A1*Xrapr;  B1Xr = B1*Xrapr;  C1Xr = C1*Xrapr;  D1Xr = D1*Xrapr;  
            A2Yr = A2*Yrapr;  B2Yr = B2*Yrapr;  C2Yr = C2*Yrapr;  D2Yr = D2*Yrapr;  
            A3Zr = A3*Zrapr;  B3Zr = B3*Zrapr;  C3Zr = C3*Zrapr;  D3Zr = D3*Zrapr;  
            tmpM = [Xrapr'*B1Xr Xrapr'*C1Xr Xrapr'*D1Xr; Yrapr'*B2Yr Yrapr'*C2Yr Yrapr'*D2Yr; Zrapr'*B3Zr Zrapr'*C3Zr Zrapr'*D3Zr];
            tmpb = [Xrapr'*A1Xr; Yrapr'*A2Yr;  Zrapr'*A3Zr];
            tilda = tmpM\tmpb;
            la = tilda(1);  ua = tilda(2); ea = tilda(3);      
        end

        ABC1 = A1-la*B1-ua*C1-ea*D1;
        ABC2 = A2-la*B2-ua*C2-ea*D2;
        ABC3 = A3-la*B3-ua*C3-ea*D3;
        r1r = ABC1*Xrapr;  % residual r1 right
        r2r = ABC2*Yrapr;  % residual r2 right
        r3r = ABC3*Zrapr;  % residual r3 right
        tmpres(step) = norm([r1r;r2r;r3r], 'fro'); 	% norm of residual
        if showinfo==1
            disp(sprintf('%3d size %4d skip %3d/%3d, conv: %3d, dist %7.1e, res: %7.1e  eig: (%10.3e%+10.3ei, %10.3e%+10.3ei, %10.3e%+10.3ei) minPr: %7.1e',...
                 step,size(U1,2)*size(U2,2)*size(U3,2),skipped,length(order),conv,dddd,tmpres(step),real(la),imag(la),real(ua),imag(ua),real(ea),imag(ea),minInnerProd))
        end
        if length(res)<step
            res(step) = tmpres(step);
        end

        if (tmpres(step)<=switcheps) && (expandrandom==0) % candidate for convergence 
            % we have a possible new eigenvalue, we refine it and test the residual again
            if forcereal
                [tilda, posmax] = max(abs(Xrapr)); Xrapr = real(Xrapr/Xrapr(posmax)); Xrapr = Xrapr/norm(Xrapr);
                [tilda, posmax] = max(abs(Yrapr)); Yrapr = real(Yrapr/Yrapr(posmax)); Yrapr = Yrapr/norm(Yrapr);
                [tilda, posmax] = max(abs(Zrapr)); Zrapr = real(Zrapr/Zrapr(posmax)); Zrapr = Zrapr/norm(Zrapr);
            end
            if refine>0
                [la,ua,ea,Xrapr2,Yrapr2,Zrapr2] = trqi_3p(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,Xrapr,Yrapr,Zrapr,refine,delta);
                RefABC1 = A1-la*B1-ua*C1-ea*D1;
                RefABC2 = A2-la*B2-ua*C2-ea*D2;
                RefABC3 = A3-la*B3-ua*C3-ea*D3;
            else
                RefABC1 = ABC1;
                RefABC2 = ABC2;
                RefABC3 = ABC3;
                Xrapr2 = Xrapr;
                Yrapr2 = Yrapr;
                Zrapr2 = Zrapr;
            end
            % residual
            r1r = RefABC1*Xrapr2;  % residual r1 right
            r2r = RefABC2*Yrapr2;  % residual r2 right
            r3r = RefABC3*Zrapr2;  % residual r3 right
            resnorm = norm([r1r;r2r;r3r], 'fro'); 	% norm of residual

            if isempty(XPr)
                primerjava = 0;
            else   
                innerprod = (Xrapr2'*B1XPl).*(Yrapr2'*C2YPl).*(Zrapr2'*D3ZPl) + (Xrapr2'*C1XPl).*(Yrapr2'*D2YPl).*(Zrapr2'*B3ZPl) ...
                          + (Xrapr2'*D1XPl).*(Yrapr2'*B2YPl).*(Zrapr2'*C3ZPl) - (Xrapr2'*B1XPl).*(Yrapr2'*D2YPl).*(Zrapr2'*C3ZPl) ...
                          - (Xrapr2'*C1XPl).*(Yrapr2'*B2YPl).*(Zrapr2'*D3ZPl) - (Xrapr2'*D1XPl).*(Yrapr2'*C2YPl).*(Zrapr2'*B3ZPl);
                primerjava = max(abs(innerprod)./InnerProds);
            end
            
            % we check the residual and inner Delta product again (as both might change because of TRQI)
            if (resnorm<delta) && (primerjava<selcrit2)
                Xrapr = Xrapr2;
                Yrapr = Yrapr2;
                Zrapr = Zrapr2;
                XPr = [XPr Xrapr];
                YPr = [YPr Yrapr];
                ZPr = [ZPr Zrapr];
                % left eigenvectors are computed for the selection criteria
                [tilda, Xlapr] = min_sing_vec(RefABC1,1);
                [tilda, Ylapr] = min_sing_vec(RefABC2,1);
                [tilda, Zlapr] = min_sing_vec(RefABC3,1);
                XPl = [XPl Xlapr];
                YPl = [YPl Ylapr];
                ZPl = [ZPl Zlapr];
                B1XPl = B1'*XPl; C1XPl = C1'*XPl; D1XPl = D1'*XPl;
                B2YPl = B2'*YPl; C2YPl = C2'*YPl; D2YPl = D2'*YPl;
                B3ZPl = B3'*ZPl; C3ZPl = C3'*ZPl; D3ZPl = D3'*ZPl;
                B1Xlapr = B1'*Xlapr; C1Xlapr = C1'*Xlapr; D1Xlapr = D1'*Xlapr;
                B2Ylapr = B2'*Ylapr; C2Ylapr = C2'*Ylapr; D2Ylapr = D2'*Ylapr;
                B3Zlapr = B3'*Zlapr; C3Zlapr = C3'*Zlapr; D3Zlapr = D3'*Zlapr;
                innerpr = (Xrapr'*B1Xlapr).*(Yrapr'*C2Ylapr).*(Zrapr'*D3Zlapr) + (Xrapr'*C1Xlapr).*(Yrapr'*D2Ylapr).*(Zrapr'*B3Zlapr) ...
                        + (Xrapr'*D1Xlapr).*(Yrapr'*B2Ylapr).*(Zrapr'*C3Zlapr) - (Xrapr'*B1Xlapr).*(Yrapr'*D2Ylapr).*(Zrapr'*C3Zlapr) ...
                        - (Xrapr'*C1Xlapr).*(Yrapr'*B2Ylapr).*(Zrapr'*D3Zlapr) - (Xrapr'*D1Xlapr).*(Yrapr'*C2Ylapr).*(Zrapr'*B3Zlapr);
                InnerProds = [InnerProds abs(innerpr)];

                minInnerProd = min(InnerProds);
                
                % new eigenvalue
                conv = conv+1;
                lambda = [lambda la];  
                mu = [mu ua];  
                eta = [eta ea];
                
                lastconv = step;
                
                noexpansion = 1;  % noexpansion, we check other Ritz vectors
                extraction = initialextraction;
                
                if showinfo>0
                    disp(sprintf('Eig (%2d): lambda: %11.4e%+11.4ei, mu: %11.4e%+11.4ei, eta: %11.4e%+11.4ei, step %4d, res: %5.1e, selcrit: %5.1e',...
                        conv,real(la),imag(la),real(ua),imag(ua),real(ea),imag(ea),step,resnorm,primerjava))
                end
            else
                if (resnorm>=delta) 
                    hist.missedTRQI = hist.missedTRQI + 1;
                    if showinfo==1
                        disp(sprintf('no TRQI convergence, step %d res0: %5.4e, res:%5.4e dist: %5.4e',step,tmpres(step),resnorm,primerjava))
                    end
                    noexpansion = 0;
                else
                    if showinfo==1
                        disp(sprintf('TRQI convergence to repeated eigenvalue (%11.4e%+11.4ei, %11.4e%+11.4ei, %11.4e%+11.4ei), step %d res0: %5.4e, res:%5.4e dist: %5.4e',...
                            real(la),imag(la),real(ua),imag(ua),real(ea),imag(ea),step,tmpres(step),resnorm,primerjava))
                    end
                    noexpansion = 1;
                    firstcand = izb + 1;
                    hist.repetead = hist.repetead  + 1;
                    selcrit1 = selcriteria(order(izb))*0.99;
                    if showinfo==1
                        disp(sprintf('selcrit lowered to %5.4e',selcrit1))
                    end
                end
            end
        else
            % no convergence in this step
            if (tmpres(step)<reschange) && (expandrandom==0) && (~strcmp(extraction,'minres'))
                if showinfo==1,  disp('Change to min. residual'),  end
                extraction='minres';
            end
            noexpansion = 0;
        end
    end % while (conv<nreig) && noexpansion
    
    % Step 2: We find new directions dx1 and dx2 for the search space 
    % --------------------------------------------------------------------------------------
    % we use orthogonal correction equation, which is 
    % preconditioned first order correction equation P2 in paper (2).                            
    if innersteps>=0
        c1 = (target(1)-la)*B1*Xrapr + (target(2)-ua)*C1*Xrapr + (target(3)-ea)*D1*Xrapr;
        c2 = (target(1)-la)*B2*Yrapr + (target(2)-ua)*C2*Yrapr + (target(3)-ea)*D2*Yrapr;
        c3 = (target(1)-la)*B3*Zrapr + (target(2)-ua)*C3*Zrapr + (target(3)-ea)*D3*Zrapr;
        
        if ~isempty(M1)
            if isa(M1,'function_handle')
                tmp1 = M1(r1r);  tmp2 = M1(c1);
            else
                tmp1 = M1*r1r;   tmp2 = M1*c1;
            end
            r1rnew = -tmp1 + tmp2*(Xrapr'*tmp1)/(Xrapr'*tmp2);
        else
            tmp2 = [];
            r1rnew = -r1r;
        end
        if ~isempty(M2)
            if isa(M2,'function_handle')
                tmp3 = M2(r2r);  tmp4 = M2(c2);
            else
                tmp3 = M2*r2r;   tmp4 = M2*c2;
            end
            r2rnew = -tmp3 + tmp4*(Yrapr'*tmp3)/(Yrapr'*tmp4);
        else
            tmp4 = [];
            r2rnew = -r2r;
        end
        if ~isempty(M3)
            if isa(M3,'function_handle')
                tmp5 = M3(r3r);  tmp6 = M3(c3);
            else
                tmp5 = M3*r3r;   tmp6 = M3*c3;
            end
            r3rnew = -tmp5 + tmp6*(Zrapr'*tmp5)/(Zrapr'*tmp6);
        else
            tmp6 = [];
            r3rnew = -r3r;
        end
        
        dxr1 = gmres_jd(tmp2, Xrapr, M1, c1, Xrapr, ABC1, Xrapr, Xrapr, r1rnew, innersteps, innertol);
        dxr2 = gmres_jd(tmp4, Yrapr, M2, c2, Yrapr, ABC2, Yrapr, Yrapr, r2rnew, innersteps, innertol);
        dxr3 = gmres_jd(tmp6, Zrapr, M3, c3, Zrapr, ABC3, Zrapr, Zrapr, r3rnew, innersteps, innertol);
        
    else
        % if innersteps<0 then we solve correction equations exactly
        dxr1 = solveCorrectionEquationExactly(ABC1,Xrapr);
        dxr2 = solveCorrectionEquationExactly(ABC2,Yrapr);
        dxr3 = solveCorrectionEquationExactly(ABC3,Zrapr);
    end
    
    % save new directions for restart
    RestX = [Xrapr+dxr1 RestX];
    RestY = [Yrapr+dxr2 RestY];
    RestZ = [Zrapr+dxr3 RestZ];
    
    % Step 3: We expand search space in directions dx1,dx2,dx3
    % --------------------------------------------------------------------------------------
    % dx1, dx2, dx3 are new directions, we need orthogonal projections on
    % orth(U1), orth(U2), and orth(U3)
    
    k = size(U1,2)+1;
    
    newxr1 = ExpandMGS(U1,dxr1);
    newxr2 = ExpandMGS(U2,dxr2);
    newxr3 = ExpandMGS(U3,dxr3);
    U1(:,k) = newxr1;
    U2(:,k) = newxr2;
    U3(:,k) = newxr3;
    
    % we compute new columns and rows of matrices
    % U1'*A1*U1,U1'*B1*U1,U1'*C1*U1,U1'*D1*U1,U2'*A2*U2,U2'*B2*U2,U2'*C2*U2, ..., U3'*D3*U3
    AU1 = [AU1 A1*newxr1];  BU1 = [BU1 B1*newxr1];  CU1 = [CU1 C1*newxr1];  DU1 = [DU1 D1*newxr1];
    AU2 = [AU2 A2*newxr2];  BU2 = [BU2 B2*newxr2];  CU2 = [CU2 C2*newxr2];  DU2 = [DU2 D2*newxr2];
    AU3 = [AU3 A3*newxr3];  BU3 = [BU3 B3*newxr3];  CU3 = [CU3 C3*newxr3];  DU3 = [DU3 D3*newxr3];

    if harmonic
       tmpx = HABC1*newxr1;
       tmpy = HABC2*newxr2;
       tmpz = HABC3*newxr3;
       HABCU1 = [HABCU1 tmpx];
       HABCU2 = [HABCU2 tmpy];
       HABCU3 = [HABCU3 tmpz];
       newh1 = ExpandMGS(Q1, tmpx);
       newh2 = ExpandMGS(Q2, tmpy);
       newh3 = ExpandMGS(Q3, tmpz);
       Q1 = [Q1 newh1];
       Q2 = [Q2 newh2];
       Q3 = [Q3 newh3];
    end
    
    step = step+1;
    
  end % while

  % Step 4: We restart J-D 
  % --------------------------------------------------------------------------------------
  % we restart with last minsize new directions ... minsize last approximations that are Delta0 distant enough distances
  if (conv<neig) && (step<maxsteps) 
     
     U1 = orth(RestX(:,1:min(size(RestX,2),minsize))); 
     U2 = orth(RestY(:,1:min(size(RestY,2),minsize))); 
     U3 = orth(RestZ(:,1:min(size(RestZ,2),minsize))); 

     if size(U1,2)<minsize, U1 = orth([U1 randn(n1,minsize-size(U1,2),class_t)]); end
     if size(U2,2)<minsize, U2 = orth([U2 randn(n2,minsize-size(U2,2),class_t)]); end
     if size(U3,2)<minsize, U3 = orth([U3 randn(n3,minsize-size(U3,2),class_t)]); end

     RestX = RestX(:,min(size(RestX,2),minsize));
     RestY = RestY(:,min(size(RestY,2),minsize));
     RestZ = RestZ(:,min(size(RestZ,2),minsize));
     
     AU1 = A1*U1;  BU1 = B1*U1;  CU1 = C1*U1;  DU1 = D1*U1;
     AU2 = A2*U2;  BU2 = B2*U2;  CU2 = C2*U2;  DU2 = D2*U2;
     AU3 = A3*U3;  BU3 = B3*U3;  CU3 = C3*U3;  DU3 = D3*U3;
     
     if harmonic
        HABCU1 = HABC1*U1;
        HABCU2 = HABC2*U2;
        HABCU3 = HABC3*U3;
        Q1 = rgs(HABCU1); 
        Q2 = rgs(HABCU2); 
        Q3 = rgs(HABCU3); 
     end
  end  
    
end % outer loop : while (conv<nreig) && (step<maxJDsteps)

% Results
if conv>0
    % we report only new eigenvectors if some were supplied in the beginning
    XR = XPr(:,end-conv+1:end); 
    YR = YPr(:,end-conv+1:end); 
    ZR = ZPr(:,end-conv+1:end); 
    XL = XPl(:,end-conv+1:end); 
    YL = YPl(:,end-conv+1:end); 
    ZL = ZPl(:,end-conv+1:end); 
else
    XR = []; YR = []; ZR = []; XL = []; YL = []; ZL = []; 
end

step = step-1;

lambda = lambda(:);
mu = mu(:);
eta = eta(:);

hist.switcheps = switcheps;
hist.minInnerProd = minInnerProd;
hist.selcrit1 = selcrit1;
hist.steps = step;
hist.InnerProds = InnerProds;

% -----------------------------------------------------------------------------------------------
% Auxiliary routine ExpandMGS
% -----------------------------------------------------------------------------------------------

function y = ExpandMGS(Q,x)

% Auxiliary routine that orthogonalizes x against the orthogonal columns 
% of U. Two steps of the repeated GS are used to maintain the stability

% Bor Plestenjak
% last revision 22.08.04

c = 0.71;
k = size(Q,2);
oldnorm = norm(x);
for j = 1:k % modified GS for additional direction
   x = x - (Q(:,j)'*x)*Q(:,j);
end
newnorm = norm(x);
if newnorm < c*oldnorm
   x = x/newnorm;
   for j = 1:k % modified GS for additional direction
      x = x - (Q(:,j)'*x)*Q(:,j);
   end
   newnorm = norm(x);
end
y = x/newnorm;

% -----------------------------------------------------------------------------------------------
% Auxiliary routine rgs (repeated Gram-Schmidt orthogonalization of A)
% -----------------------------------------------------------------------------------------------

function Q = rgs(A)

c=0.71;

[m,n] = size(A);

for i = 1:n
   Q(:,i) = A(:,i);
   oldnorm = norm(Q(:,i));
   r = oldnorm;
   for j = 1:i-1
      r = Q(:,j)'*Q(:,i);
      Q(:,i) = Q(:,i)-r*Q(:,j);
   end
   newnorm = norm(Q(:,i));
   if newnorm < c*oldnorm 
       for j = 1:i-1
          r = Q(:,j)'*Q(:,i);
          Q(:,i) = Q(:,i)-r*Q(:,j);
       end
       newnorm = norm(Q(:,i));
   end
   Q(:,i) = Q(:,i)/newnorm;
end

% -----------------------------------------------------------------------------------------------
% Auxiliary routine gmres_jd
% -----------------------------------------------------------------------------------------------

function [x, relres, k] = gmres_jd(A1x, A1y, A2, A3x, A3y, A4, A5x, A5y, b, maxit, tol)

%GMRES_JD is a version of gmres that is used in correction equations 
%   in Jacobi-Davidson type methods for two-parameter eigenvalue problems 
%   (can be applied to other similar JD methods as well)
%   The matrix vector multiplicaton in each step is 
%   y = P1* A2 * P3 * A4 * P5 * x,
%   where 
%      P1 : a) P1 is projection (I-inv(A1y'*A1x)*A1x*A1y') or 
%           b) P1=I when A1x=[] 
%      A2 : a) A2=A2 or 
%           b) A2=I when A2A=[] or
%           c) A2 is given by function_handle that evaluates A2(x)
%      P3 : a) P3 is projection (I-inv(A3y'*A3x)*A3x*A3y') or 
%           b) P3=I when A3x=[] 
%      A4 : a) A4=A4 or 
%           b) A4=I when A4A=[] or
%      P5 : a) P5 is projection (I-inv(A5y'*A5x)*A5x*A5y') or 
%           b) P5=I when A5x=[] 
%
%   function [x, relres, k] = gmres_jd(A1x, A1y, A2, A3x, A3y, A4, A5x, A5y, b, maxit, tol)
%   Assumes x0 = 0
%   In:
%      b     : righthand side
%      maxit : maximum number of steps (no restarts)
%      tol   : tolerance
%   Out:
%      x     : gmres solution
%      relres: obtained residual reduction
%      k     : number of steps performed

% Last revision: 31.08.2015
% Adapted from gmres_fast and mgs code by Michiel Hochstenbach
% Bor Plestenjak

class_t = superiorfloat(A1x,A1y,A2,A3x,A3y,A4,A5x,A5y,b);

c = 0.71; % ~= 1/sqrt(2)
H = zeros(maxit+1, maxit, class_t);
V = zeros(length(b), maxit, class_t);

if (maxit < 1) || (tol >= 1)
   relres = 1;
   x = b;
   return
end 

if ~isempty(A1x),  A1d = inv(A1y'*A1x);  else A1d = 1;  end
if ~isempty(A3x),  A3d = inv(A3y'*A3x);  else A3d = 1;  end
if ~isempty(A5x),  A5d = inv(A5y'*A5x);  else A5d = 1;  end

rho0 = norm(b);
b = b/rho0;

v = b;
Gamma = 1;
k = 0;
rho = 1;
if tol == 0
	tol0 = Inf;
else
	tol0 = 1/(tol*tol); 
end

while (k < maxit) && (rho < tol0) 
   k = k+1;
   V(:,k) = v;
  
   % multiplication by A, projections, and preconditioner
   % ----------------------------------------------------
   if ~isempty(A5x),  v = v - A5x*(A5d*(A5y'*v));  end
   if ~isempty(A4),   v = A4*v;  end
   if ~isempty(A3x),  v = v - A3x*(A3d*(A3y'*v));  end
   if ~isempty(A2)
      if isa(A2,'function_handle')
          v = A2(v);
      else
          v = A2*v;
      end
   end
   if ~isempty(A1x),  v = v - A1x*(A1d*(A1y'*v));  end
   % ----------------------------------------------------
 
   % Arnoldi step 
 
   H(k+1,k) = 1;
   l = 1;
   norm_input = norm(v);
   norm_output = 0;
   while (l <= 2) && (norm_output < c * norm_input)
      for j = 1:k
         inpr = V(:,j)'*v;
         v = v - inpr*V(:,j);    
         H(j,k) = H(j,k) + H(k+1,k) * inpr;
      end
	  norm_output = norm(v);
      v = v/norm_output;
      H(k+1,k) = H(k+1,k) * norm_output;
      l = l+1;
   end
   
   gamma = H(k+1,k);
  
   if gamma == 0 % Lucky break-down
      break
   else
      gamma = -Gamma*H(1:k,k)/gamma; 
      Gamma = [Gamma gamma];
      rho   = rho + gamma'*gamma;
   end
end
if gamma == 0   % Lucky break-down
   relres = 0;
   y = zeros(k,1, class_t);
   y(1) = 1;
   x = V(:,1:k)*(H(1:k,1:k) \ y);
else            % solve in least square sense 
   relres = 1 / sqrt(rho);
   y = zeros(k+1,1, class_t);
   y(1) = 1;
   x = V(:,1:k)*(H(1:k+1,1:k) \ y);
end

x = rho0 * x;

% -----------------------------------------------------------------------------------------------
% Auxiliary routine solveCorrectionEquationExactly
% -----------------------------------------------------------------------------------------------

function y = solveCorrectionEquationExactly(A,u)

% exact solution of correction equation (I-u*u')*A*(I-u*u')y = -r, 
% where r = A u 
%
% Here A = A0 - theta I, such that A0 u - theta u is orthogonal to u

% matrix might be close to singular, we turn off the warning
warning off
z = A\u;
warning on
a = 1/(z'*u);
y = -u + a*z;
