function [labels] = dbscanclust(reach , coreDist, ordering, eps)
    %dbscanclust - Cluster the reachability of OPTICS with the dbscan
    %cutoff method

    labels = zeros(size(reach,1),1);
    
    far = reach > eps;
    prox = coreDist <= eps;
    
    labels(ordering) = cumsum( far(ordering) & prox(ordering))-1;
    labels(far & ~prox) = -1;
    
end
