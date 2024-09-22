clear
close all
clc

%% load matrices 
fileContent = load('matrices_brassTeflon.mat');
E0 = fileContent.E0;
E1 = fileContent.E1;
E2 = fileContent.E2;
M  = fileContent.M;
R  = fileContent.R;
c  = fileContent.c;

typeCoupling = 'S';                                                         % type of coupling (one solid)
attThreshold = 4000;                                                        % maximum attenuation to consider
w = 2*pi*linspace(1e-2, 7, 281).';                                          % frequency

%% solver
% allocate arrays of horizontal and vertical wavenumbers
% arrays are allocated as complex nan to avoid extra zeros in plots 
k = nan(length(w), 2^3*size(E0,1))*(1+1i);                                  % horizontal wavenumber
kyL = k;                                                                    % vertical wavenumber, bottom, pressure wave
kyS = k;                                                                    % vertical wavenumber, bottom, pressure wave
for i = 1:length(w)
    kappa = w(i)./c;
    [~, allEV] = eig_Leaky_all(E0,E1,-E2,M,R,typeCoupling,kappa,w(i),[]);
    nSol = size(allEV,1);
    k(i,1:nSol)   = allEV(:,1);
    kyL(i,1:nSol) = allEV(:,2);
    kyS(i,1:nSol) = allEV(:,3);
end

att = imag(k)*20/log(10)*1000;                                              % attenuation

%% filter 
% filter out incoming waves and negative attenuation
indRemove = (real(kyL)<0) | (real(kyS)<0) | (att>attThreshold) | (att<0);
k(indRemove) = nan + 1i*nan;
att(indRemove) = nan;

%% plot
load reference_BrassTeflon.mat
plotResults(w,k,att,refCp,refAtt,attThreshold);





