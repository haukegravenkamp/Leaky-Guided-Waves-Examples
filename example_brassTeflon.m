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

typeCoupling = 'S';             % type of coupling (one solid)
attThreshold = 4000;            % maximum attenuation to consier
w = 2*pi*linspace(1e-3, 7, 281).'; % frequency

%% solver
n = size(E0,1);
NN = nan(length(w), 2^3*n);
k = NN + 1i*NN;
kyL = k;
kyS = k;
for i = 1:length(w)
    kappa = w(i)./c;
    [lambda, tmp_lambda] = eig_Leaky_all(E0,E1,-E2,M,R,typeCoupling,kappa,w(i),[]);
    k(i,1:numel(lambda))  = lambda;
    kyL(i,1:numel(lambda))  = tmp_lambda(:,2);
    kyS(i,1:numel(lambda))  = tmp_lambda(:,3);
end

att = imag(k)*20/log(10)*1000;   % attenuation

%% filter 
indRemove = (real(kyL)<0) | (real(kyS)<0) | (att>attThreshold) | (att<0);
k(indRemove) = nan + 1i*nan;
att(indRemove) = nan;

%% plot
load reference_BrassTeflon.mat
plotResults(w,k,att,refCp,refAtt,attThreshold);




