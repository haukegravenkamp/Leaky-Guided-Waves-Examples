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

typeCoupling = 'FS';                % type of coupling (one solid)
attThreshold = 2000;                % maximum attenuation to consier
w = 2*pi*linspace(1e-3, 3, 121).';  % frequency

%% solver
n = size(E0,1);
NN = nan(length(w), 2^3*n);
k = NN + 1i*NN;
kyBp = k;
kyBs = k;
kyT = k;
for i = 1:length(w)
    kappa = w(i)./c;
    [lambda, tmp_lambda] = eig_Leaky_all(E0,E1,-E2,M,R,typeCoupling,kappa,w(i),[]);
    k(i,1:numel(lambda))  = lambda;
    kyBp(i,1:numel(lambda))  = tmp_lambda(:,2);
    kyBs(i,1:numel(lambda))  = tmp_lambda(:,3);
    kyT(i,1:numel(lambda))  = tmp_lambda(:,4);
end

att = imag(k)*20/log(10)*1000;   % attenuation

%% filter
indRemove = (real(kyBp)>-1e-2) | (real(kyBs)>-1e-2) | (real(kyT)<1e-2) | (att>attThreshold) | (att<1e-2);
k(indRemove) = nan + 1i*nan;
att(indRemove) = nan;

%% plot
load reference_compositeTeflonOil.mat
plotResults(w,k,att,refCp,refAtt,attThreshold)






