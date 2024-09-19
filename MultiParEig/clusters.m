function [order,clstart,clsize,newx] = clusters(x,tol)

%CLUSTERS  Groups elements into clusters
%
% [order,clstart,clsize,newx] = CLUSTERS(x,eps)
% orders elements of x into clusters based on their distance
%
% x(i) and x(j) are in the same cluster if abs(x(i)-x(j))/(1+abs(x(i))) <= tol
% Output:
%   order : permutation of the initial x that reorders x into clusters
%   clstart : indices of starts of the clusters
%   clsize : sizes of clusters
%   newx : reodered x where each cluster is represented by the average

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% First revision: 8.9.2015
% BP 11.11.2023: support for infinite values

n = size(x,1);
% change infinite values into finite
isfx = isfinite(x);
fn = sum(isfx);
xx = x;
if fn==0
    % all values are infinite and in the same cluster
    order = 1:n;
    clstart = 1;
    clsize = n;
    newx = x;
    return
elseif fn<n
    bigx = 10*sum(abs(x(isfinite(x))));
    xx(not(isfx)) = bigx;
end
dif = abs((xx*ones(1,n,class(x))-ones(n,1,class(x))*xx.')*(diag(1./(1+abs(xx))))) < tol;
order = [];
clstart = [];
clsize = [];
free = 1:n;
while any(free)
   cluster = [];
   [tilda1,tilda2,elem] = find(free); %#ok<*ASGLU>
   newindices = elem(1);
   while ~isempty(newindices)
	  free(newindices) = 0;
   	  cluster = [cluster newindices]; %#ok<*AGROW>
      [p,q] = find(dif(cluster,:));
      q = q(:)';
      d = q((1:end-1)')==q((2:end)');  
      q(d) = [];   
      [tlda1,tilda2,tmpfree] = find(free);
      newindices = sort([q tmpfree]);
      tmp3 = newindices(1:end-1)==newindices(2:end);
      newindices = newindices(tmp3);
   end
   clstart = [clstart length(order)+1];
   order = [order cluster];
   clsize = [clsize length(cluster)];
end   

newx = x(order,1);
for k = [clstart; clsize]
   if isfinite(newx(clstart))
      newx(k(1):(k(1)+k(2)-1),1) = sum(newx(k(1):(k(1)+k(2)-1),1))/k(2);
   end
end


