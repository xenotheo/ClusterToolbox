function [labels] = spectralClustering(X, N, varargin)
% 
% Spectral Clustering  
%         The function takes the N dimmensional points and
%         the number of desired clusters as inputs and
%         returns the label for each point. This algorithm
%         will try to find the N smallest eigenvalues of 
%         Similarity matrix using eigs(...). In case of a 
%         badly conditioned Similarity matrix, it will use
%         svds(...) instead and the results may be inaccurate.
% SYNTAX
%
%   [ L ] = spectralClustering( X, N )
%   [ L ] = spectralClustering( X , N, 'h', 2)
%
% INPUT
% 
%   X       The dataset of N dimmensional points                [N-by-d]
%           Or A Sparse Graph of Pairwise distances             [N-by-N]
%
%   N       Number of desired clusters                          [Integer]
% OPTIONAL
%   H       Bandwidth parameter of similar Matrix               [Scalar]
%           {default:1}
% OUTPUT
% 
%   LABELS  An array contaning the labels for each point.      [N-vector]
%           
% EXAMPLE
%   load clustDatasets.mat
%   X  = clustdataset{7}.Dataset;
%   L  = spectralClustering( X, 15, 'h', 0.5);
%   nC = length(unique(L));
%   figure
%   scatter(X(:,1),X(:,2), [], L);
%   colormap(hsv(nC));
%   title('Clustering Result of Spectral Clustering')



%% PARSE OPTIONAL INPUT
opt = parseOptArgs( varargin{:} );

%% CALCULATE LAPLACIAN MATRIX
if issparse(X)
    W =  spfun(@(x) exp( - (x/opt.h).^2 ),X);
    W(1:size(W,1)+1 :end) = 0;
    Diag = diag(sum(W,2));
    L = (Diag.^0.5)\W/(Diag.^0.5);
else
    W = pdist2(X,X);
    W = exp( - (W./(opt.h)).^2 );
    W(1:size(W,1)+1 :end) = 0;
    Diag = diag(sum(W,2));
    L = (Diag.^0.5)\W/(Diag.^0.5); 
end

%% CALCULATE K LARGEST EIGENVECTORS
V = calculateEigen(L,N);



%% CLUSTER EXTRACTION WITH K-MEANS
labels = kmeans(V , N);




%% LOCAL FUNCTION: CALCULATE EIGENVECTORS

function V = calculateEigen(L,N)
try 
    [V , D] = svds( L, N,'largest');
    D = diag(D);
    [~ , I]= sort(D,'descend');
    V = V(:,I);
    V = bsxfun(@(x,y) x./y ,V,sqrt(sum(V.^2,2)));

catch ME 
    if strcmp(ME.identifier,'MATLAB:eigs:SingularA')
        fprintf('\n...The Similarity matrix is badly conditioned\n');
        fprintf('Using SVDS instead. Results may be inaccurate\n');
        [V , D] = svds( L, N, 'largest');
        D = diag(D);
        [~ , I]= sort(D,'descend');
        V = V(:,I);
        V = bsxfun(@(x,y) x./y ,V,sqrt(sum(V.^2,2)));
    else
        rethrow(ME)
    end
end 
    
    
end

%% LOCAL FUNCTION: OPTIONAL ARGUMENT PARSING
%
function opt = parseOptArgs(varargin)

% default values
dflt.h             = 0.25;
%dflt.method        = 'U';
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