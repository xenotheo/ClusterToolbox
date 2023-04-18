function [rho, delta, minIdx, perm] = genDecisionGraph(D, dc, type)
% GENDECISIONGRAPH - Generate decision graph, given input points
%   
% SYNTAX
%
%   [RHO, DELTA, IDX, PERM] = GENDECISIONGRAPH( DIST, DC )
%
% INPUT
%
%   DIST
%   DC
%
% OUTPUT
%
%   RHO
%   DELTA
%   IDX
%   PERM
%
% DESCRIPTION
%
%   GENDECISIONGRAPH generates the decision graph, as explained in [1].
%
% REFERENCES
%   
%   [1] Rodriguez, A. and Laio, A., 2014. Clustering by fast search
%       and find of density peaks. Science, 344(6191), pp.1492-1496.
%
% See also      k-means
%

  % get input size
  n = size( D, 1 );
  
  % make sure same point is not accounted for
  if issparse(D)
    D = D + spdiags( inf( n, 1 ),0,n,n );
    
    switch lower( type )
    
    case {'hard', 'uniform', 'box'}
      rho = sum( spfun(@(x) x < dc,D), 2 );
      
    case {'soft', 'gaussian'}
      rho = sum( spfun(@(x) exp( - (x/dc).^2 ),D), 2 );
      
    otherwise
      error( 'Unknown type' )
      
    end
  else
    D = D + diag( inf( n, 1 ) );
    switch lower( type )
    
    case {'hard', 'uniform', 'box'}
      rho = sum( D < dc, 2 );
      
    case {'soft', 'gaussian'}
      rho = sum( exp( - (D/dc).^2 ), 2 );
      
    otherwise
      error( 'Unknown type' )
      
    end
  end
  
  
  % select type of kernel
  
      
  % distance to higher distance (for each point)
  [~, perm] = sort(rho, 'descend');
  
  % maximum distance in feature space
  % permutation of the pair-wise distances (by density, descend)
  D = D(perm, perm);
  
  % remove distances to larger points
  D = tril( D );
  
  
  if(issparse(D))
  

     
      [ii,~] = find(D);
      delta = accumarray(ii,nonzeros(D),[],@min);
      minIdx = zeros(n,1);
      for i =1:n
          minIdx(i) = find(D(i,:)==delta(i),1,'first');
      end
      delta(1) = max(D(spfun(@(x) ~isinf(x),D)));
      delta(isinf(delta))=delta(1);
  else
      % do not account for zero distances
      D( D <= 0 ) = Inf;
      [delta,minIdx] = min(D,[],2);
       
      % change Inf of first distance
      delta(1) = max(D(~isinf(D)));
      
      % permute back
      
      
  end
      delta( perm ) = delta;
      minIdx( perm) = perm( minIdx );
end


%%------------------------------------------------------------
%
%
% AUTHORS
%
%   Dimitris Floros                         fcdimitr@auth.gr
%  
% VERSION
%   0.3 - March 15, 2023
%   
%
% CHANGELOG
%   0.3 (March 15, 2023) - Xenofon
%       *added sparse functionality 
%
%
%   0.2 (Jun 07, 2018) - Dimitris
%       * updated distance of max density point
% 
%   0.1 (Apr 25, 2018) - Dimitris
%       * initial implementation
%
% ------------------------------------------------------------

