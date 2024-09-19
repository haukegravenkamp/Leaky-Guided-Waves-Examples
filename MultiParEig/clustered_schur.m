function [Q,S,order,start,csize,lambda] = clustered_schur(A,epscluster)

%CLUSTERED_SCHUR  Schur decomposition with clustering 
%
% [Q,S,order,start,csize,lambda] = CLUSTERED_SCHUR(A,epscluster)

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 23.10.2023: initial version

[tmpQ,tmpT] = schur(A,'complex');
tmplambda = diag(tmpT);
[order,start,csize,lambda] = clusters(tmplambda,epscluster);
n = length(order);
reverse = order(n:-1:1);
invord = zeros(1,n);
for j = 1:length(reverse)
    invord(reverse(j)) = j;
end
[Q,S] = ordschur(tmpQ,tmpT,invord);