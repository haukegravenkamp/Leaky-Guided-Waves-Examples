function [lambda,mu,XR,YR,XL,YL,res,hist] = twopareigs_jd(A1,B1,C1,A2,B2,C2,neig,opts)

%TWOPAREIGS_JD   Jacobi-Davidson method for a two-parameter eigenvalue problem
%
% [lambda,mu,XR,YR,XL,YL,res] = TWOPAREIGS_JD(A1,B1,C1,A2,B2,C2,neig,opts)
% returns neig eigenvalues of the two-parameter eigenvalue problem
%
% A1*x = lambda*B1*x + mu*C1*x
% A2*y = lambda*B2*y + mu*C2*y
%
% using the Jacobi-Davidson method.
%
% Input:
%   - A1, B1, C1, A2, B2, C2 : matrices
%   - neig : number of eigenvalues (6)             
%   - opts : options (see below)
%
% Output:
%   - lambda, mu : eigenvalue parts (eigenvalues are (lambda(j),mu(j))
%   - XR, YR : components of decomposable right eigenvectors (eigenvector is kron(XR(:,j),YR(:,j)), such that
%       (A1-lambda(j)*B1-mu(j)*C1)*XR(:,j)=0
%       (A2-lambda(j)*B2-mu(j)*C2)*YR(:,j)=0
%   - XL, YL : components of decomposable left eigenvectors (eigenvector is kron(XL(:,j),YL(:,j)), such that
%       (A1-lambda(j)*B1-mu(j)*C1)'*XL(:,j)=0
%       (A2-lambda(j)*B2-mu(j)*C2)'*YL(:,j)=0
%   - res : norms of residuals for all J-D steps, try plot(log10(res))
%   - hist : additional information for developers
%
% Options are (default values in parenthesis):
%   - delta : convergence criterion (norm of the residual) for the outer iteration (1e4*eps*max inf norm of the matrices)
%   - switcheps : criterion when (norm of the residual) to switch to TRQI refinement for the new eigenpair candidate (1e2*delta)
%   - X0r,Y0r : initial search space (rand(n1,1) and rand(n2,1))
%   - maxsize : maximum dimension of search spaces before restart (10)
%   - minsize : dimension of search spaces after restart (5)
%   - maxsteps : maximum number of outer iteration steps (100)
%   - innersteps : number of GMRES steps for the correction equation (5, set to -1 to solve the correction equation exactly)
%   - innertol : tolerance in GMRES for the correction equation (0)
%   - window : number of Ritz values with the smallest |mu| that we compute, set to 0 (default) for all Ritz values
%   - extraction : extraction method, choices are: 
%        'maxdist' : the maximum distance from the target,
%        'mindist' : the minimum distance from the target (default),
%        'minres'  : the smallest residual,
%        'maxlambda': the eigenvalue with lambda with the largest real part
%        'minmu' : the eigenvalue with the smallest |mu|
%   - showinfo : display nothing (default), 1: all values, 2: just eigenvalues
%   - target : target for the eigenvalues ([0 0])
%   - reschange : switch to minimum residual extraction when residual norm is small (0) (no effect if smaller than switcheps)
%   - rcsteps : switch to minimum residual if no eigenvalues are found in rcsteps iteration (Inf)
%   - XPr,XPl,YPl,YPr  - prohibited directions (previously computed left and right eigenvectors) - ([],[],[],[])
%   - M1, M2 : left preconditioners - inverses for A1,B1,C1 and A2,B2,C2, respectively, can also be a function_handle that returns M1*x1 or M2*x2
%   - harmonic : set to 1 to use harmonic instead of Ritz values (0) - use this for interior eigenvalue
%   - forcereal : set to 1 if you know that eigenvectors and eigenvalues are real, default is 0
%   - refine : (3) number of TRQI steps to refine an eigenpair candidate
%   - solveropts : options for twopareigs of twopareig ([])
%   - fp_type: numeric type to use ('single', 'double', or 'mp' (needs MCT),
%     use only if you want to change the default - the superior type of input data
%
% See also: THREEPAREIGS_JD, TWOPAREIG, TWOPAREIGS, TWOPAREIGS_IRA, TWOPAREIGS_KS,
% TWOPAREIGS_SI

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
%      Subspace methods for 3-parameter eigenvalue problems, arXiv:1802:07386

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 06.11.2016 : Switch to TRQI refinement and some small changes
% BP Fixed bug in the preconditioned first order correction 
% BP 31.01.208 : new features from JD for 3-parameter problems : exact solution, trqi refinement, new restart
% Modified for MCT

% Last revision: 31.01.2018

narginchk(6,8)

if nargin<7, neig = 6; end
if nargin<8, opts = []; end
if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(A1,B1,C1,A2,B2,C2);
end

% Make sure all inputs are of the same numeric type.
if ~isa(A1,class_t), A1 = numeric_t(A1,class_t); end;
if ~isa(B1,class_t), B1 = numeric_t(B1,class_t); end;
if ~isa(C1,class_t), C1 = numeric_t(C1,class_t); end;
if ~isa(A2,class_t), A2 = numeric_t(A2,class_t); end;
if ~isa(B2,class_t), B2 = numeric_t(B2,class_t); end;
if ~isa(C2,class_t), C2 = numeric_t(C2,class_t); end;

n1 = size(A1,1); 
n2 = size(A2,1);

EmptyMatrix = numeric_t([],class_t);

if isfield(opts,'delta')
    delta = opts.delta;            
else
    maxnor = max([norm(A1,'inf') norm(B1,'inf') norm(C1,'inf') norm(A2,'inf') norm(B2,'inf') norm(C2,'inf')]);
    delta = numeric_t('1e4*eps',class_t)*maxnor;            
end
if isfield(opts,'switcheps'),   switcheps = opts.switcheps;    else switcheps = numeric_t('1e2',class_t)*delta;   end
if isfield(opts,'minsize'),     minsize = opts.minsize;        else minsize = 5;             end
if isfield(opts,'maxsize'),     maxsize = opts.maxsize;        else maxsize = 10;            end
if isfield(opts,'maxsteps'),    maxsteps = opts.maxsteps;      else maxsteps = 100;          end
if isfield(opts,'showinfo'),    showinfo = opts.showinfo;      else showinfo = 2;            end
if isfield(opts,'target'),      target = opts.target;          else target = [0 0];          end
if isfield(opts,'X0r'),         X0r = opts.X0r;                else X0r = rand(n1,1,class_t); end
if isfield(opts,'Y0r'),         Y0r = opts.Y0r;                else Y0r = rand(n2,1,class_t); end
if isfield(opts,'innersteps'),  innersteps = opts.innersteps;  else innersteps = 5;          end
if isfield(opts,'innertol'),    innertol = opts.innertol;      else innertol = 0;            end
if isfield(opts,'window'),      window = opts.window;          else window = 0;              end
if isfield(opts,'extraction'),  extraction = opts.extraction;  else extraction = 'mindist';  end
if isfield(opts,'reschange'),   reschange = opts.reschange;    else reschange = 0;           end % 10^(-2.5);   end
if isfield(opts,'rcsteps'),     rcsteps = opts.rcsteps;        else rcsteps = Inf;           end
if isfield(opts,'XPr'),         XPr = opts.XPr;                else XPr = EmptyMatrix;       end
if isfield(opts,'XPl'),         XPl = opts.XPl;                else XPl = EmptyMatrix;       end
if isfield(opts,'YPr'),         YPr = opts.YPr;                else YPr = EmptyMatrix;       end
if isfield(opts,'YPl'),         YPl = opts.YPl;                else YPl = EmptyMatrix;       end
if isfield(opts,'M1'),          M1 = opts.M1;                  else M1 = [];                 end
if isfield(opts,'M2'),          M2 = opts.M2;                  else M2 = [];                 end
if isfield(opts,'harmonic'),    harmonic = opts.harmonic;      else harmonic = 0;            end
if isfield(opts,'forcereal'),   forcereal = opts.forcereal;    else forcereal = 0;           end
if isfield(opts,'refine'),      refine = opts.refine;          else refine = 3;              end
if isfield(opts,'solveropts'),  solveropts = opts.solveropts;  else solveropts = [];         end
if isfield(opts,'selcrit1'),    selcrit1 = opts.selcrit1;      else selcrit1 = 1e-1;         end
if isfield(opts,'selcrit2'),    selcrit2 = opts.selcrit2;      else selcrit2 = 1e-4;         end

if ~isfield(solveropts,'usesparse'),  solveropts.usesparse = 0; end;

if ~isempty(XPl)
   B1XPl = B1'*XPl; C1XPl = C1'*XPl;
   B2YPl = B2'*YPl; C2YPl = C2'*YPl;
else
   B1XPl = EmptyMatrix; C1XPl = EmptyMatrix; 
   B2YPl = EmptyMatrix; C2YPl = EmptyMatrix; 
end

% Initial search spaces and other initialization
U1 = rgs(X0r); 
U2 = rgs(Y0r); 
RestX = EmptyMatrix; 
RestY = EmptyMatrix;

initialextraction = extraction;
lambda = EmptyMatrix; 
mu = EmptyMatrix; 
conv = 0; % no convergence yet
step = 1; % number of steps
maxsteps = maxsteps + 1; % as we start counting with 1 

lastconv = 0; % step with the last convergence

AU1 = A1*U1;  BU1 = B1*U1;  CU1 = C1*U1; 
AU2 = A2*U2;  BU2 = B2*U2;  CU2 = C2*U2; 

if harmonic
    HABC1 = A1-target(1)*B1-target(2)*C1;
    HABC2 = A2-target(1)*B2-target(2)*C2;
    HABCU1 = HABC1*U1;
    HABCU2 = HABC2*U2;
    Q1 = rgs(HABCU1); 
    Q2 = rgs(HABCU2); 
end    

% selection criteria for Delta0 orthogonality
InnerProds = EmptyMatrix;
if ~isempty(XPr) 
   for j = 1:size(XPr,2)
        innerpr = (XPr(:,j)'*B1XPl(:,j)).*(YPr(:,j)'*C2YPl(:,j)) - (XPr(:,j)'*C1XPl(:,j)).*(YPr(:,j)'*B2YPl(:,j));
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
        SA1 = Q1'*HABCU1; SB1 = Q1'*BU1; SC1 = Q1'*CU1;
        SA2 = Q2'*HABCU2; SB2 = Q2'*BU2; SC2 = Q2'*CU2;
    else
        SA1 = U1'*AU1;    SB1 = U1'*BU1; SC1 = U1'*CU1;
        SA2 = U2'*AU2;    SB2 = U2'*BU2; SC2 = U2'*CU2;
    end
    
    if window>0
        [Al,Au,AXr,AYr] = twopareigs(SA1,SB1,SC1,SA2,SB2,SC2,window,solveropts);
    else
        [Al,Au,AXr,AYr] = twopareig(SA1,SB1,SC1,SA2,SB2,SC2,solveropts);
    end
    
    tmpvel = length(Al);

    if harmonic
        Al = Al + target(1);
        Au = Au + target(2);
    end

    noexpansion = 1;   
    firstcand = 1;
    
    % in a loop we can extract more than one eigenpair from the subspace
    while (conv<neig) && noexpansion 
        
       switch extraction
           case 'maxdist'  % Ritz value farthest to the target
               ritznorm = abs(Al-target(1)).^2 + abs(Au-target(2)).^2;
               [tilda,order] = sort(-ritznorm);
           case 'minmu'  % Ritz value with the smallest |mu|
               [tilda,order] = sort(abs(Au));
           case 'mindist'  % Ritz value closest to the target 
               ritznorm = abs(Al-target(1)).^2 + abs(Au-target(2)).^2;
               [tilda,order] = sort(ritznorm);
           case 'maxlambda'  % Ritz value with the largest real part
               [tilda,order] = sort(-Al);
           case 'minres'  % Ritz pair with the smallest residual               
               rn = zeros(length(Al),1);
               for kk = 1:length(Al)
                   ABCU1 = AU1-Al(kk)*BU1-Au(kk)*CU1;
                   ABCU2 = AU2-Al(kk)*BU2-Au(kk)*CU2;
                   r1r = ABCU1*AXr(:,kk); % residual r1 right
                   r2r = ABCU2*AYr(:,kk); % residual r2 right
                   rn(kk,1) = norm([r1r;r2r], 'fro'); % norm of residual
               end   
               [tilda,order] = sort(rn);
           otherwise
               error('Unknown extraction option')
        end
         
        % we select the Ritz vector that is Delta0 orthogonal to computed eigenpairs
        if isempty(XPr)
            selcriteria = zeros(length(Al),1);
        else   
            SAXr = U1*AXr;  SAYr = U2*AYr;
            innerprod = (SAXr'*B1XPl).*(SAYr'*C2YPl) - (SAXr'*C1XPl).*(SAYr'*B2YPl);
            primerjava = abs(innerprod)./(ones(tmpvel,1)*InnerProds);
            selcriteria = max(primerjava,[],2);
        end
        
        izb = firstcand;
        while selcriteria(order(izb))>selcrit1
            izb = izb+1;
            if izb>tmpvel  break;  end
        end
        skipped = izb-1;
   
        if izb>tmpvel  
            % all directions are not enough "Delta0 orthogonal", so we pick the one with the minimal distance    
            [tilda,pos] = min(selcriteria);
            izb = 1;
            expandrandom = 1;  % safety flag that prevents taking already computed eigenvalue
        else
            pos = order(izb);
            expandrandom = 0;
        end;
        
        dddd = selcriteria(pos); 
        la = Al(pos);  ua = Au(pos);
        lapr(step) = la;  uapr(step) = ua;  % best Ritz values 
        Xrapr = U1*AXr(:,pos);  Yrapr = U2*AYr(:,pos);  % best Ritz vectors
        
        if harmonic
            % TRQ improvement
            A1Xr = A1*Xrapr;  B1Xr = B1*Xrapr;  C1Xr = C1*Xrapr; 
            A2Yr = A2*Yrapr;  B2Yr = B2*Yrapr;  C2Yr = C2*Yrapr; 
            tmpM = [Xrapr'*B1Xr Xrapr'*C1Xr; Yrapr'*B2Yr Yrapr'*C2Yr];
            tmpb = [Xrapr'*A1Xr; Yrapr'*A2Yr];
            tilda = tmpM\tmpb;
            la = tilda(1);  ua = tilda(2);
            lapr(step) = la;  uapr(step) = ua; % corrected Ritz values 
        end

        ABC1 = A1-la*B1-ua*C1;
        ABC2 = A2-la*B2-ua*C2;
        r1r = ABC1*Xrapr;  % residual r1 right
        r2r = ABC2*Yrapr;  % residual r2 right
        tmpres(step) = norm([r1r;r2r], 'fro'); 	% norm of residual
        if showinfo==1
            disp(sprintf('%3d size %4d skip %3d/%3d, conv: %3d, dist %7.1e, res: %7.1e  eig: (%10.3e%+10.3ei,%10.3e%+10.3ei) minPr: %7.1e',...
                 step,size(U1,2)*size(U2,2),skipped,tmpvel,conv,dddd,tmpres(step),real(la),imag(la),real(ua),imag(ua),minInnerProd))
        end
        if length(res)<step
            res(step) = tmpres(step);
        end

        if (tmpres(step)<=switcheps) && (expandrandom==0) % candidate for convergence 
            % we have a possible new eigenvalue, we refine it and test the residual again
            if forcereal
                [tilda, posmax] = max(abs(Xrapr)); Xrapr = real(Xrapr/Xrapr(posmax)); Xrapr = Xrapr/norm(Xrapr);
                [tilda, posmax] = max(abs(Yrapr)); Yrapr = real(Yrapr/Yrapr(posmax)); Yrapr = Yrapr/norm(Yrapr);
            end
            if refine>0
                [la,ua,Xrapr,Yrapr] = trqi(A1,B1,C1,A2,B2,C2,Xrapr,Yrapr,refine,delta);
                ABC1 = A1-la*B1-ua*C1;
                ABC2 = A2-la*B2-ua*C2;
            end
            % residual
            r1r = ABC1*Xrapr;  % residual r1 right
            r2r = ABC2*Yrapr;  % residual r2 right
            resnorm = norm([r1r;r2r], 'fro'); 	% norm of residual

            if isempty(XPr)
                primerjava = 0;
            else   
                innerprod = (Xrapr'*B1XPl).*(Yrapr'*C2YPl) - (Xrapr'*C1XPl).*(Yrapr'*B2YPl);
                kvocienti = abs(innerprod)./InnerProds;
                primerjava = max(abs(innerprod)./InnerProds);
            end
            
            % we check the residual and inner Delta product again (as both might change because of TRQI)
            if (resnorm<delta) && (primerjava<selcrit2)
                XPr = [XPr Xrapr];
                YPr = [YPr Yrapr];
                % left eigenvectors are computed for the selection criteria
                [tilda, Xlapr] = min_sing_vec(ABC1,1);
                [tilda, Ylapr] = min_sing_vec(ABC2,1);
                XPl = [XPl Xlapr];
                YPl = [YPl Ylapr];
                B1XPl = B1'*XPl; C1XPl = C1'*XPl;
                B2YPl = B2'*YPl; C2YPl = C2'*YPl;
                B1Xlapr = B1'*Xlapr; C1Xlapr = C1'*Xlapr;
                B2Ylapr = B2'*Ylapr; C2Ylapr = C2'*Ylapr;
                innerpr = (Xrapr'*B1Xlapr).*(Yrapr'*C2Ylapr) - (Xrapr'*C1Xlapr).*(Yrapr'*B2Ylapr);
                InnerProds = [InnerProds abs(innerpr)];

                minInnerProd = min(InnerProds);
                
                % new eigenvalue
                conv = conv+1;
                lambda = [lambda la];  
                mu = [mu ua];  
                extraction = initialextraction;
                
                lastconv = step;
                if izb < tmpvel
                    firstcand = izb + 1;
                    noexpansion = 1;  % noexpansion, we check other Ritz vectors
                else
                    noexpansion = 0;
                end
                
                if showinfo>0
                    disp(sprintf('Eig (%2d): lambda: %+11.4e%+11.4ei, mu:%11.4e%+11.4ei, step %4d, res: %7.1e, selcrit: %7.1e',...
                        conv,real(la),imag(la),real(ua),imag(ua),step,resnorm,primerjava))
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
                        disp(sprintf('TRQI convergence to repeated eigenvalue (%+5.4e %+5.4ei, %+5.4e %+5.4ei), step %d res0: %5.4e, res:%5.4e dist0: %5.4e dist: %5.4e',...
                            real(la),imag(la),real(ua),imag(ua),step,tmpres(step),resnorm,selcriteria(order(izb)),primerjava))
                    end
                    noexpansion = 1;
                    firstcand = izb + 1;
                    hist.repetead = hist.repetead  + 1;
                    selcrit1 = selcriteria(order(izb))*0.99;
                    if showinfo==1
                        disp(sprintf('selcrit lowered to %5.4e',selcrit1))
                    end
                end
                if (hist.missedTRQI>1) && (mod(hist.missedTRQI,5)==1) && (switcheps>500*delta)
                    switcheps = switcheps/2;
                    if showinfo==1
                        disp(sprintf('Aggressive delta decreased to %5.4e, step %d',switcheps,step))
                    end
                end
            end
        else
            % no convergence in this step
            if ((tmpres(step)<reschange) || (step-lastconv>rcsteps)) && (expandrandom==0) && (~strcmp(extraction,'minres'))
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
        c1 = (target(1)-la)*B1*Xrapr + (target(2)-ua)*C1*Xrapr;
        c2 = (target(1)-la)*B2*Yrapr + (target(2)-ua)*C2*Yrapr;
        
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
        dxr1 = gmres_jd(tmp2, Xrapr, M1, c1, Xrapr, ABC1, Xrapr, Xrapr, r1rnew, innersteps, innertol,class_t);
        dxr2 = gmres_jd(tmp4, Yrapr, M2, c2, Yrapr, ABC2, Yrapr, Yrapr, r2rnew, innersteps, innertol,class_t);
        
    else
        % if innersteps<0 then we solve correction equations exactly
        dxr1 = solveCorrectionEquationExactly(ABC1,Xrapr);
        dxr2 = solveCorrectionEquationExactly(ABC2,Yrapr);
    end
    
    % save new directions for restart
    RestX = [Xrapr+dxr1 RestX];
    RestY = [Yrapr+dxr2 RestY];
    
    % Step 3: We expand search space in directions dx1,dx2,dx3
    % --------------------------------------------------------------------------------------
    % dx1 and dx2 are new directions, we need orthogonal projections on orth(U1) and orth(U2)
    
    k = size(U1,2)+1;
    
    newxr1 = ExpandMGS(U1,dxr1);
    newxr2 = ExpandMGS(U2,dxr2);
    U1(:,k) = newxr1;
    U2(:,k) = newxr2;
    
    % we compute new columns and rows of matrices U1'*A1*U1,U1'*B1*U1,...,U2'*C2*U2
    AU1 = [AU1 A1*newxr1];  BU1 = [BU1 B1*newxr1];  CU1 = [CU1 C1*newxr1]; 
    AU2 = [AU2 A2*newxr2];  BU2 = [BU2 B2*newxr2];  CU2 = [CU2 C2*newxr2]; 

    if harmonic
       tmpx = HABC1*newxr1;
       tmpy = HABC2*newxr2;
       HABCU1 = [HABCU1 tmpx];
       HABCU2 = [HABCU2 tmpy];
       newh1 = ExpandMGS(Q1, tmpx);
       newh2 = ExpandMGS(Q2, tmpy);
       Q1 = [Q1 newh1];
       Q2 = [Q2 newh2];
    end
    
    step = step+1;
    
  end % while

  % Step 4: We restart J-D 
  % --------------------------------------------------------------------------------------
  % we restart with last minsize new directions ... minsize last approximations that are Delta0 distant enough distances
  if (conv<neig) && (step<maxsteps) 
     
     U1 = orth(RestX(:,1:min(size(RestX,2),minsize))); 
     U2 = orth(RestY(:,1:min(size(RestY,2),minsize))); 

     if size(U1,2)<minsize, U1 = orth([U1 randn(n1,minsize-size(U1,2),class_t)]); end
     if size(U2,2)<minsize, U2 = orth([U2 randn(n2,minsize-size(U2,2),class_t)]); end

     RestX = RestX(:,min(size(RestX,2),minsize));
     RestY = RestY(:,min(size(RestY,2),minsize));
     
     AU1 = A1*U1; BU1 = B1*U1; CU1 = C1*U1; 
     AU2 = A2*U2; BU2 = B2*U2; CU2 = C2*U2; 
     
     if harmonic
        HABCU1 = HABC1*U1;
        HABCU2 = HABC2*U2;
        Q1 = rgs(HABCU1); 
        Q2 = rgs(HABCU2); 
     end
  end  
    
end % outer loop : while (conv<nreig) && (step<maxJDsteps)

% Results
if conv>0
    % we report only new eigenvectors if some were supplied in the beginning
    XR = XPr(:,end-conv+1:end); 
    YR = YPr(:,end-conv+1:end); 
    XL = XPl(:,end-conv+1:end); 
    YL = YPl(:,end-conv+1:end); 
else
    XR = []; YR = []; XL = []; YL = [];
end    
step = step-1;

lambda = lambda(:);
mu = mu(:);

hist.switcheps = switcheps;
hist.minInnerProd = minInnerProd;
hist.selcrit1 = selcrit1;
hist.steps = step;

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

function [x, relres, k] = gmres_jd(A1x, A1y, A2, A3x, A3y, A4, A5x, A5y, b, maxit, tol, class_t)

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
