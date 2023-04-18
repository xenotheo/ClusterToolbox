function [params] = setparams(dataset , algorithm)

  dist = pdist2(dataset.Dataset , dataset.Dataset);  
  nClust = length(unique(dataset.Cluster));
  meanD  = mean(mean(dist));
  switch  algorithm 
    
    case {'kmeans'}
        params = {dataset.Dataset , nClust};
            
    case {'kmedoids'}
        params = {dataset.Dataset , nClust,...
                 'Algorithm','Clara'};    
    case {'affinityPropagation'}
        params = {dataset.Dataset, 'preference', -max(max(dist))};
    case {'cluster'}
        params = {'MaxClust',nClust};
        
    case {'dbscan'}
        params = {dataset.Dataset,'eps', calculateOptDist(dataset.Dataset,4) , 'minPts' , 4};
    
    case {'optics'}
        optDist = calculateOptDist(dataset.Dataset,4);
        params = {dataset.Dataset , 4 ,'epsilon', 1.425*optDist, 'epsdbscan', optDist};
    case {'spectralClustering'}  
        params = {dataset.Dataset , nClust, 'h', 1};
    case {'meanshift'} 
        params = {dataset.Dataset, (sqrt(2)/3)*sqrt(meanD)};
    case {'densityCluster'} 
        params = {dataset.Dataset, (sqrt(2)/3)*sqrt(meanD), 'nClusters', nClust};
    
    
        
  end





end