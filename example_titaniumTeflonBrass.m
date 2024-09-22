clear
close all
clc

%% load matrices
fileContent = load('matrices_titaniumTeflonBrass.mat');
E0 = fileContent.E0;
E1 = fileContent.E1;
E2 = fileContent.E2;
M  = fileContent.M;
R  = fileContent.R;
c  = fileContent.c;

typeCoupling = 'SS';                                                        % type of coupling (solid-solid)
attThreshold = 30000;                                                       % maximum attenuation to consider
w = 2*pi*linspace(1e-3, 10, 201).';                                         % frequency

%% solver
% allocate arrays of horizontal and vertical wavenumbers
% arrays are allocated as complex nan to avoid extra zeros in plots 
k = nan(length(w), 2^3*size(E0,1))*(1+1i);                                  % horizontal wavenumber
kyBp = k;                                                                   % vertical wavenumber, bottom, pressure wave
kyBs = k;                                                                   % vertical wavenumber, bottom, shear wave
kyTp = k;                                                                   % vertical wavenumber, top, pressure wave
kyTs = k;                                                                   % vertical wavenumber, top, shear wave
for i = 1:length(w)
    kappa = w(i)./c;
    [~, allEV] = eig_Leaky_all(E0,E1,-E2,M,R,typeCoupling,kappa,w(i),[]);
    nSol = size(allEV,1);
    k(i,1:nSol)     = allEV(:,1);
    kyBp(i,1:nSol)  = allEV(:,2);
    kyBs(i,1:nSol)  = allEV(:,3);
    kyTp(i,1:nSol)  = allEV(:,4);
    kyTs(i,1:nSol)  = allEV(:,5);
end

att = imag(k)*20/log(10)*1000;                                              % attenuation

%% filter
% filter out incoming waves and negative attenuation
indRemove = (real(kyBp)>-1e-2) | (real(kyBs)>-1e-2) | (real(kyTp)<1e-2)...
    | (real(kyTs)<1e-2) | (att>attThreshold) | (att<1e-2);
k(indRemove) = nan + 1i*nan;
att(indRemove) = nan;

%% plot
load reference_titaniumTeflonBrass.mat
plotResults(w,k,att,refCp,refAtt,attThreshold)






