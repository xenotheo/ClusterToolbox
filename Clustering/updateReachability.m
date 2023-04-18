function [reachability , predecessor ] = updateReachability(coreDist, reachability, predecessor,...
                                 point_index,    processed, D ,  nbrs)
%OPTICS Helper Function
%Calculates and Updates the reachability distance of Points

   unprocessed = nbrs(~processed(nbrs));

   if (isempty(unprocessed)) 
       return
   end
   [~,~,D] = find(D);
   reachnew = max(full(D) , coreDist * ones(1,length(D)));
   perm = bsxfun...
              (@lt,...
                  reachnew,reachability(unprocessed)');
                        
   reachability(unprocessed(perm)) = reachnew(perm);
   predecessor(unprocessed(perm))  = point_index;
end