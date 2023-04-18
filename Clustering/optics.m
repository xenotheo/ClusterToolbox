function [labels, reachability, order, coreDist, predecessor] = optics(X , minPts , varargin)
% 
% OPTICS  Matlab implementation
%         The function takes the N dimmensional points as inputs and
%         returns the label for each point.
%   
% SYNTAX
%
%   [ L, R, O ] = optics( X,  NPTS )
%   [ L, R, O ] = optics( X , NPTS, 'epsilon', 3, 'epsdbscan', 2)
%
% INPUT
% 
%   X          The dataset of N dimmensional points                [N-by-d]
%
%   MINPTS     Minimum number of points required                  [Integer]
%              for a point to be considered as core. 
%          
% OPTIONAL
%   EPS        The neighborhood radius value                       [Scalar]
%              {default:5}
%
%   EPSDBSCAN  The radius cutoff value used by dbscan clustering   [Scalar]
%              on reachability. 
%              {default:1.5}
%
%   VIZ        Set it to true and the function reachabilityPlot   [Boolean]
%              
%              will be evoked. It will plot the reachability/ordering 
%              diagram as well as the clustering result of dbscan
%              cut off method. I highly suggest to set this to true.
%              {default:false}
%             
%   
% OUTPUT
% 
%   LABELS       An array contaning the labels for each point.     [N-by-1]
%                -1 means noise point.
%   
%   REACHABILITY An array contaning the reachability distance      [N-by-1]
%                for each point. Seed points have Inf reachability
%
%   ORDER        An array contaning the sequence that the points   [N-by-1]
%                were prossed.
%
%   COREDIST     An array contaning the distance of the            [N-by-1]
%                'minPts'th NN. If a point doesn't have minPts
%                within a radius of EPS it have Inf coreDist       
%
%   PREDECESSOR  In i_th-row contains the index of the point       [N-by-1]
%                that called the i_th point
% EXAMPLE
%   load clustDatasets.mat
%   X          = clustdataset{2}.Dataset;
%   [L, R, D]  = optics(X,8 ,'epsilon',2,'viz',true,'epsdbscan',1.25);
  

%% PARSE OPTIONAL INPUT
opt = parseOptArgs( varargin{:} );


%% RANGE SEARCH

[idxCol,dist] = rangesearch(X, X, opt.epsilon);
dist = cellfun(@(x) x+eps,dist,'UniformOutput',false);
nNbr = cellfun( @(x) numel(x), idxCol );
n = numel( idxCol );

% row indices (for sparse matrix formation convenience)
idxRow = arrayfun( @(n,i) i * ones( 1, n ), nNbr, (1:n)', ...
                   'UniformOutput', false );
  
% sparse matrix formation
D = sparse( [idxRow{:}], [idxCol{:}], [dist{:}], n, n );

%% 1: CORE POINTS

idxCore = sum(D>0, 2) >= minPts;
coreDist = Inf * ones(n,1);

coreDist(idxCore) = D( sub2ind ([n,n], find(idxCore) , ...
                    cellfun(@(x) x(minPts),idxCol(idxCore,:))));

             
reachability = Inf*ones(n,1);
predecessor = int64(-1*ones(n,1));
processed = zeros(n, 1);
order = uint64(zeros(n,1));

%% 2: CALCULATE REACHABILITY , ORDERING
for i=1:n
  pivot = find(processed ==0);
  
  [~ , I] = min(reachability(pivot));
  point = pivot(I);
  processed(point) = true;
  order(i) = point;
  
  if (coreDist(point)~=Inf)
      
      [reachability, predecessor] = updateReachability(coreDist(point), reachability, predecessor,...
                         point, processed,D(point,~processed),  idxCol{point});
   
  end
end

%% 3. EXTRACT LABELS VIA DBSCAN CUTOFF
labels = dbscanclust(reachability , coreDist, order, opt.epsdbscan);

if (opt.viz)
    reachabilityPlot(reachability , order, opt.epsdbscan, X, labels);
end

%% LOCAL FUNCTION: OPTIONAL ARGUMENT PARSING
%
function opt = parseOptArgs(varargin)

    % default values
    dflt.epsilon       = 5;
    dflt.epsdbscan     = 1.6;
    dflt.viz           = false;

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