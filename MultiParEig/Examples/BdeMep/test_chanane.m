% demo boucherif & chanane

% computes values from B. Chanane & A. Boucherif : Two-parameter
% Sturm-Liouville problems

a = 0;
b = 0.7;
c = 1;
p = -1;
q = 0;
r = 0;
s = 1;
t = @(x) x;
bc = [1 0; 1 0];

n1 = 40;
n2 = 40;
opts = [];

[z1,A1,B1,C1,G1,k1,r1] = bde2mep(a,b,p,q,r,s,t,bc,n1,opts);
[z2,A2,B2,C2,G2,k2,r2] = bde2mep(b,c,p,q,r,s,t,bc,n1,opts);

[lambda,mu] = twopareigs(A1,B1,C1,A2,B2,C2,100);
ind = intersect(find(lambda>0),find(mu>0));
mu1 = sqrt(lambda(ind));
mu2 = sqrt(mu(ind));

sols = real([mu1 mu2])


 