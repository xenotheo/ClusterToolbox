function [clustCent, data2cluster] = meanShiftClust(y, h, scale)


  if ~exist('scale', 'var') || isempty( scale )
    scale = 1;
  end

  % get number of points
  n = size( y, 1 );

  % one cluster is there for sure
  clustCent    = y(1,:);
  data2cluster = [1 zeros( 1, n-1 ) ];

  % identify more clusters -- scan over the query points
  for i = 2:n
    
    % get current point
    y_temp = y(i,:);
    
    % find distance from cluster centers
    d = pdist2( y_temp, clustCent );
    
    % find clusters nearby
    closeCenters = d < h/2 .* scale;
    
    if any( closeCenters )          % if a nearby cluster exists
      
      % add this point to the cluster
      clustIdx = find( closeCenters, 1, 'first' );
      
      % update the cluster center
      clustCent(clustIdx,:) = mean( [clustCent(clustIdx,:); y_temp], 1 );
      
    else                            % otherwise, create new
      
      % create new cluster and add this point
      clustCent(end+1, :) = y_temp; %#ok
      clustIdx = size( clustCent, 1 );
      
    end
    
    % update cluster for point i
    data2cluster(i) = clustIdx;
    
  end

end