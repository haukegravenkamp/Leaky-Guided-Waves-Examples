clear
close all
clc

%% load matrices
fileContent = load('matrices_brassWater.mat');
E0 = fileContent.E0;
E1 = fileContent.E1;
E2 = fileContent.E2;
M  = fileContent.M;
R  = fileContent.R;
c  = fileContent.c;

typeCoupling = 'F';                                                         % type of coupling (one fluid)
attThreshold = 2000;                                                        % maximum attenuation to consider
w = 2*pi*linspace(0.001, 4, 300).';                                         % frequency

%% solver
% allocate arrays of horizontal and vertical wavenumbers
k = nan(length(w), 2^3*size(E0,1))*(1+1i);                                  % horizontal wavenumber
kyB = k;                                                                    % vertical wavenumber, bottom, pressure wave
for i = 1:length(w)
    kappa = w(i)./c;
    [~, allEV] = eig_Leaky_all(E0,E1,-E2,M,R,typeCoupling,kappa,w(i),[]);
    nSol = size(allEV,1);
    k(i,1:nSol)   = allEV(:,1);
    kyB(i,1:nSol) = allEV(:,2);
end
kyT = -kyB;

att = imag(k)*20/log(10)*1000;                                              % attenuation

%% filter
indRemove = (real(kyB)>1e-2) | (real(kyT)<-1e-2) | (att>attThreshold) | (att<-1e-2);
k(indRemove) = nan + 1i*nan;
att(indRemove) = nan;

%% plot
load reference_BrassWater.mat
[figC, figA] = plotResults(w,k,att,refCp,refAtt,attThreshold,{},{'LineStyle','none','Marker','.'});