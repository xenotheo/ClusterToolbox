function halo = findHalo(D, rho, cid, d_c)
% FINDHALO - Find cluster halo (noise)
%   
% SYNTAX
%
%   [ HALO ] = FINDHALO( DIST, RHO, CID, DC )
%
% INPUT
%
%   DIST
%   RHO
%   CID
%   DC
%
% OUTPUT
%
%   HALO
%
% DESCRIPTION
%
%   FINDHALO finds the points of the dataset that belong to the halo.
%
% DEPENDENCIES
%
%   <none>
%
%
% See also      findClusterID, genDecisionGraph
%
  
  n = length( rho );
  k = length( unique( cid ) );

  % find for each cluster the border density
  halo = zeros( n, 1 );
  
  for i = 1:k

    % all density of cluster i
    m_i = (cid == i);
    d = rho(m_i);

    % all dist within d_c (cluster i -- other clusters)
    m_j = max( (D(m_i, ~m_i) <= d_c ), [], 2);

    % density threshold (maximum density of )
    thres = max( d( m_j ) );
    if (~isempty(thres))
      halo( m_i ) = thres;
    end

  end

  % filter out all the points < den_thres
  halo = rho < halo;
  
end


%%------------------------------------------------------------
%
% AUTHORS
%
%   Dimitris Floros                         fcdimitr@auth.gr
%
% VERSION
%
%   0.1 - April 27, 2018
%
% CHANGELOG
%
%   0.1 (Apr 27, 2018) - Dimitris
%       * initial implementation
%
% ------------------------------------------------------------

