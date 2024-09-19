function [r,choice]= numrank(d, opts)

%NUMRANK   Find a gap in singular values
%
% [r,choice] = NUMRANK(d, opts) returns rank from a set of singular values d
%
% Arguments: 
%   d : a set of singular values (e.g., result of the function svd)
%   opts: options
% 
% Options in opts:
%   - rankeps: 1e-12 (discard all singular values smaller from rankeps*initialnorm or rankeps*d(1))
%   - initialnorm: 0 (if nonzero then rank is the sum of values larger than rankeps*initialnorm)
%   - showrank: 0 (display diagnostics)
%   - fixedrank: -1 (if nonnegative, rank is prescribed)
%
% Output:
%   r : numerical rank
%   choice : how was rank selected (0: larger than rankeps*d(1), -1: fixedrank)

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% P. Holoborodko, Advanpix LLC.
% FreeBSD License, see LICENSE.txt

% BP 16.09.2015: new options fixedrank and fixedrankgap, report choice, mingap -> mingapdif
% PH 22.11.2016: precision-independent version
% BP 26.11.2016: minor change in the printout
% PH 26.11.2016: code simplifications and clean-ups.
% BP 28.11.2023: simplified version without heuristics, we compare singular values to the 
%                norms of the initial Delta matrices, the only exception is prescribed rank

class_t = class(d);

if nargin<2, opts=[]; end
if isfield(opts,'rankeps'),      rankeps      = opts.rankeps;       else rankeps      = numeric_t('1e-12',class_t);  end
if isfield(opts,'showrank'),     showrank     = opts.showrank;      else showrank     = 0;                           end
if isfield(opts,'fixedrank'),    fixedrank    = opts.fixedrank;     else fixedrank    = -1;                          end
if isfield(opts,'initialnorm'),  initialnorm =  opts.initialnorm;   else initialnorm  = 0;                           end

% These options might be set by the calling function for additional info
% that is displayed if showrank > 0 (the names of the calling function and original matrix size)
if isfield(opts,'call'),         call = opts.call;                  else call = '';        end
if isfield(opts,'m1'),           m1 = opts.m1;                      else m1 = 0;           end
if isfield(opts,'m2'),           m2 = opts.m2;                      else m2 = 0;           end

n = length(d);
choice = 0;

if n == 0
    % size 0, rank is 0
    r = 0; 
    choice = 0;
    if showrank
        fprintf('(%5d x %5d) r: %5d  gap: %8.1e d(1): %5.1e d(r-1): %5.1e d(r): %5.1e d(r+1): %5.1e d(n): %5.1e ver: %d %s rankeps: %5.1e\n',...
          m1,m2,r,0,0,0,0,0,0,0,call,rankeps);
    end
    return 
end

if initialnorm>0
    rankeps = rankeps*initialnorm;
else
    rankeps = rankeps*d(1);
end

% fixed rank
if fixedrank>=0
    r = fixedrank;
    choice = -1;
else
    r = sum(d > rankeps);
    choice = 0;
end
    
if showrank
    if r == 0
        dprev = 0;
        dran = 0;
        dnext = d(1);
        gap = log10(d(1)/rankeps);
    elseif r == n
        dprev = d(n-1);
        dran = d(n);
        dnext = 0;
        gap = log10(d(n)/rankeps);
    elseif r == 1
        dprev = 0;
        dran = d(r);
        if r<n
            dnext = d(r+1);
            gap = log10(d(1)/d(2));
        else
            dnext = 0;
            gap = log10(d(1)/rankeps);
        end
    else
        dprev = d(r-1);
        dran = d(r);
        dnext = d(r+1);
        gap = log10(d(r)/d(r+1));
    end
    fprintf('(%5d x %5d) r: %5d  gap: %8.1e d(1): %5.1e d(r-1): %5.1e d(r): %5.1e d(r+1): %5.1e d(n): %5.1e ver: %d %s rankeps: %5.1e\n',...
        m1,m2,r,gap,d(1),dprev,dran,dnext,d(n),choice,call,rankeps);
end


