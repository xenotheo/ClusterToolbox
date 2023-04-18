function [labels , c, y] = meanshift(x, h, varargin)
% MEANSHIFT - Mean shift implementation
%   
% SYNTAX
%
%   [L]        = MEANSHIFT( XIN, BAND )
%   [L, C , Y] = MEANSHIFT( ..., 'epsilon', EPSILON )
%   [L, C , Y] = MEANSHIFT( ..., 'verbose', VERBOSE )
%   [L, C , Y] = MEANSHIFT( ..., 'display', DISPLAY )
%
% INPUT
%
%   XIN         Input data (for clustering)                 [n-by-d]
%   BAND        Bandwidth value                             [scalar]
%   
% OPTIONAL
% 
%   EPSILON     Threshold for convergence                   [scalar]
%               {default: 1e-4*h}
%   VERBOSE     Print iteration number & error?             [boolean]
%               {default: false}
%   DISPLAY     Plot results of each iteration?             [boolean]
%               (only for 2D points)
%               {default: false}
%
% OUTPUT
%   
%   L           Label for each point after doing a revision  [n-by-1]
%               on the shifted Values Y.
%  
%   C           Co-ordinates of cluster centers              [C-by-d]   
%
%   Y           Final points location after mean shift       [n-by-d]
%
% DESCRIPTION
%
%   YOUT = MEANSHIFT(XIN,BAND) implements mean shift algorithm on
%   input points XIN, using Gaussian kernel with bandwidth BAND.
%   The local maxima of each point is then recorded in the output
%   array YOUT.
%
% DEPENDENCIES
%
%   <none>
%
% LOCAL-FUNCTIONS
% 
%   rangesearch2sparse
%   parseOptArgs
%
% See also      kmeans
%
  
  %% PARAMETERS

  % stoping threshold
  opt.epsilon = 1e-4*h;
  opt.verbose = false;
  opt.display = false;
  
  
  %% PARSE OPTIONAL INPUTS
  
  opt = parseOptArgs(opt, varargin{:});
  
  
  %% INITIALIZATION
  
  % number of points -- dimensionality
  [n, d] = size( x );
  
  % initialize output points to input points
  y = x;
  
  % mean shift vectors (initialize to infinite)
  m = inf;
  
  % iteration counter
  iter = 0;
  
  if opt.display && d == 2
    fig = figure(1337);
    set(fig, 'name', 'real_time_quiver')
  end
  
  while norm(m,'fro') > opt.epsilon  % --- iterate unitl convergence
  
    iter = iter + 1;
    
    % find pairwise distance matrix (inside radius)
    [I, D] = rangesearch( x, y, h );
    D      = cellfun( @(x) x.^2, D, 'UniformOutput', false );
    W      = rangesearch2sparse( I, D );
    
    % compute kernel matrix
    W = spfun( @(x) exp( -x / (2*h^2) ), W );

    % make sure diagonal elements are 1
    W = W + spdiags( ones(n,1), 0, n, n );
    
    % compute new y vector
    y_new = W * x;
    
    % normalize vector
    y_new = y_new ./ sum( W, 2 );
        
    % calculate mean-shift vector
    m = y_new - y;
    
    if opt.display && d == 2
      
      figure(1337)
      clf
      hold on
      scatter( y(:,1), y(:,2) );
      quiver( y(:,1), y(:,2), m(:,1), m(:,2), 0 );
      pause(0.3)
      
    end
    
    % update y
    y = y_new;
    
    if opt.verbose
      fprintf( ' Iteration %d - error %.2g\n', iter, norm(m,'fro') );
    end    
    
  end % while (m > epsilon)
  
    [c , labels] = meanShiftClust(y,h);
end
  


%% LOCAL FUNCTION: CREATE SPARSE MATRIX FROM RANGE SEARCH

function mat = rangesearch2sparse(idxCol, dist)
% INPUT         idxCol  Index columns for matrix        [n-cell]
%               dist    Distances of points             [n-cell]
% OUTPUT        mat     Sparse matrix with distances    [n-by-n sparse]

  % number of neighbors for each point
  nNbr = cellfun( @(x) numel(x), idxCol );
    
  % number of points
  n = numel( idxCol );

  % row indices (for sparse matrix formation convenience)
  idxRow = arrayfun( @(n,i) i * ones( 1, n ), nNbr, (1:n)', ...
                     'UniformOutput', false );
  
  % sparse matrix formation
  mat = sparse( [idxRow{:}], [idxCol{:}], [dist{:}], n, n );

  
end



%% LOCAL FUNCTION: PARSE OPTIONAL ARGUMENTS

function opt = parseOptArgs (dflt, varargin)
% INPUT         dflt    Struct with default parameters  [struct]
%               <name-value pairs>                      [varargin]
% OUTPUT        opt     Updated parameters              [struct]
  
  %% INITIALIZATION
    
  ip = inputParser;
    
  ip.CaseSensitive   = false;
  ip.KeepUnmatched   = false;
  ip.PartialMatching = true;
  ip.StructExpand    = true;
  
  
  %% PARAMETERS
  
  argNames = fieldnames( dflt );
  for i = 1 : length(argNames)
    addParameter( ip, argNames{i}, dflt.(argNames{i}) );
  end
  
  
  %% PARSE AND RETURN
  
  parse( ip, varargin{:} );
  
  opt = ip.Results;
  
  
  %% SET EMPTY VALUES TO DEFAULTS
  
  for i = 1 : length(argNames)
    if isempty( opt.(argNames{i}) )
      opt.(argNames{i}) = dflt.(argNames{i});
    end
  end
  
end