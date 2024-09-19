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

typeCoupling = 'SS';                 % type of coupling (one solid)
attThreshold = 30000;                % maximum attenuation to consier
w = 2*pi*linspace(1e-3, 10, 201).';  % frequency

%% solver
n = size(E0,1);
NN = nan(length(w), 2^3*n);
k = NN + 1i*NN;
kyBp = k;
kyBs = k;
kyTp = k;
for i = 1:length(w)
    kappa = w(i)./c;
    [lambda, tmp_lambda] = eig_Leaky_all(E0,E1,-E2,M,R,typeCoupling,kappa,w(i),[]);
    k(i,1:numel(lambda))  = lambda;
    kyBp(i,1:numel(lambda))  = tmp_lambda(:,2);
    kyBs(i,1:numel(lambda))  = tmp_lambda(:,3);
    kyTp(i,1:numel(lambda))  = tmp_lambda(:,4);
    kyTs(i,1:numel(lambda))  = tmp_lambda(:,5);
end

att = imag(k)*20/log(10)*1000;   % attenuation

%% filter
indRemove = (real(kyBp)>-1e-2) | (real(kyBs)>-1e-2) | (real(kyTp)<1e-2)  | (real(kyTs)<1e-2) | (att>attThreshold) | (att<1e-2);
k(indRemove) = nan + 1i*nan;
att(indRemove) = nan;

%% plot
load reference_titaniumTeflonBrass.mat
plotResults(w,k,att,refCp,refAtt,attThreshold)






