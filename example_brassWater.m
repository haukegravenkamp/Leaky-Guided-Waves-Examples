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

typeCoupling = 'F';                 % type of coupling (one solid)
attThreshold = 2000;                % maximum attenuation to consier
w = 2*pi*linspace(0.001, 4, 300).'; % frequency

%% solver
n = size(E0,1);
NN = nan(length(w), 2^3*n);
k = NN + 1i*NN;
kyB = k;
kyS = k;
for i = 1:length(w)
    kappa = w(i)./c;
    [lambda, tmp_lambda] = eig_Leaky_all(E0,E1,-E2,M,R,typeCoupling,kappa,w(i),[]);
    k(i,1:numel(lambda))  = lambda;
    kyB(i,1:numel(lambda))  = tmp_lambda(:,2);
end
kyT = -kyB;

att = imag(k)*20/log(10)*1000;   % attenuation

%% filter
indRemove = (real(kyB)>1e-2) | (real(kyT)<-1e-2) | (att>attThreshold) | (att<-1e-2);
k(indRemove) = nan + 1i*nan;
att(indRemove) = nan;

%% plot
load reference_BrassWater.mat
[figC, figA] = plotResults(w,k,att,refCp,refAtt,attThreshold,{},{'LineStyle','none','Marker','.'});






