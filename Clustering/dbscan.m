function [labels] = dbscan(X, varargin)
% 
% dbScan  matlab implementation
%         The function takes the N dimmensional points as inputs and
%         returns the label for each point.
%   
% SYNTAX
%
%   [ L ] = dbscan( X )
%   [ L ] = dbscan( X ,'eps', 2, 'minPts', 10)
%
% INPUT
% 
%   X       The dataset of N dimmensional points                [N-by-d]
%   
% OPTIONAL
%   EPS     The neighborhood radius value                       [Scalar]
%           {default:2}
%
%   MINPTS  Minimum number of points required                   [Scalar]
%           to form a dense region. 
%           {default:10}
%
% OUTPUT
% 
%   LABELS  An array contaning the labels for each point.      [N-vector]
%           0 means noise point.
% EXAMPLE
%   load clustdataset.mat
%   X  = clustdataset{7}.Dataset;
%   L  = dbscan(X,'eps', 0.4, 'minPts', 10)
%   nC = length(unique(L));
%   figure
%   scatter(X:,1),X(:,2), [], L)
%   colormap(hsv(nC))
%   hold on
%   scatter(X(L==0,1),X(L==0,2), [], 'k')
%   title('Clustering Result of DBSCAN')
%   hold off

%% PARSE OPTIONAL INPUT
opt = parseOptArgs( varargin{:} );



%% RANGE SEARCH

[idxCol,dist] = rangesearch(X, X, opt.eps);

nNbr = cellfun( @(x) numel(x), idxCol );
n = numel( idxCol );

% row indices (for sparse matrix formation convenience)
idxRow = arrayfun( @(n,i) i * ones( 1, n ), nNbr, (1:n)', ...
                   'UniformOutput', false );
  
% sparse matrix formation
W = sparse( [idxRow{:}], [idxCol{:}], [dist{:}], n, n );

%% 1: CORE POINTS

idxCore = sum(W>0, 2) >= opt.minPts;

%% 2: COMPONENTS

bins = conncomp(graph(W(idxCore, idxCore)));

%% 3: OTHERS

idxOthers = ~idxCore;

[i,j,v] = find(W(idxOthers, idxCore));
S = accumarray( [i bins(j)'], 1./v, ...
                [full(sum(idxOthers)), max(bins)], ...
                @max, [], false );

[val,cidOthers] = max(S, [], 2);

cidOthers(val == 0) = 0;

%% WRAP-UP

% cluster membership vector
labels = zeros(n, 1);

% set core points and others
labels(idxCore)   = bins;
labels(idxOthers) = cidOthers;



%% LOCAL FUNCTION: OPTIONAL ARGUMENT PARSING
%
function opt = parseOptArgs(varargin)

% default values
dflt.eps           = 2;
dflt.minPts        = 10;
dflt.distPreComp   = false;
dflt.metric        = 'euclidean';

% inputParser initialization
ip = inputParser;

ip.CaseSensitive   = false;
ip.KeepUnmatched   = false;
ip.PartialMatching = true;
ip.StructExpand    = true;

% set parameters & default values
argNames = fieldnames( dflt );
for ii = 1 : length(argNames)
  addParameter( ip, argNames{ii}, dflt.(argNames{ii}) );
end

% parse input
parse( ip, varargin{:} );
opt = ip.Results;

% use default values for any options set to []
for ii = 1 : length(argNames)
  if isempty( opt.(argNames{ii}) )
    opt.(argNames{ii}) = dflt.(argNames{ii});
  end
end
end
end