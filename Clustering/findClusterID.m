function cid = findClusterID(idxMin,idxC,ord)
% FINDCLUSTERID - Find cluster memberships
%   
% SYNTAX
%
%   CID = FINDCLUSTERID( IDXMIN, IDXC, ORD )
%
% INPUT
%
%   IDXMIN
%   IDXC
%   N
%
% OUTPUT
%
%   CID
%
% DESCRIPTION
%
%   FINDCLUSTERID assigns each point with the corresponding cluster
%   center.
%
% DEPENDENCIES
%
%   
%
%
% See also      
%
  
  n = length( idxMin );
  k = length( idxC );
  
  cid = zeros( n, 1 );
  
  cid( idxC ) = 1:k;
  
  for i = ord(:)'
    if (isempty(find(idxC == i, 1)))
      cid(i) = cid( idxMin(i) );
    end
  end
  
  
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

