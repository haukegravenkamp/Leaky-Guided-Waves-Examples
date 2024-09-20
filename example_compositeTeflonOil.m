clear
close all
clc

%% load matrices
fileContent = load('matrices_compositeTeflonOil.mat');
E0 = fileContent.E0;
E1 = fileContent.E1;
E2 = fileContent.E2;
M  = fileContent.M;
R  = fileContent.R;
c  = fileContent.c;

typeCoupling = 'FS';                                                        % type of coupling (fluid/solid)
attThreshold = 2000;                                                        % maximum attenuation to consider
w = 2*pi*linspace(1e-3, 3, 121).';                                          % frequency

%% solver
% allocate arrays of horizontal and vertical wavenumbers
k = nan(length(w), 2^3*size(E0,1))*(1+1i);                                  % horizontal wavenumber
kyBp = k;                                                                   % vertical wavenumber, bottom, pressure wave
kyBs = k;                                                                   % vertical wavenumber, bottom, shear wave
kyTp = k;                                                                   % vertical wavenumber, top, pressure wave
for i = 1:length(w)
    kappa = w(i)./c;
    [~, allEV] = eig_Leaky_all(E0,E1,-E2,M,R,typeCoupling,kappa,w(i),[]);
    nSol = size(allEV,1);
    k(i,1:nSol)    = allEV(:,1);
    kyBp(i,1:nSol) = allEV(:,2);
    kyBs(i,1:nSol) = allEV(:,3);
    kyTp(i,1:nSol) = allEV(:,4);
end

att = imag(k)*20/log(10)*1000;                                              % attenuation

%% filter
indRemove = (real(kyBp)>-1e-2) | (real(kyBs)>-1e-2) | (real(kyTp)<1e-2) | (att>attThreshold) | (att<1e-2);
k(indRemove) = nan + 1i*nan;
att(indRemove) = nan;

%% plot
load reference_compositeTeflonOil.mat
plotResults(w,k,att,refCp,refAtt,attThreshold)






