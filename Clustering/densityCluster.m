function [ clusterID, idxCenter, halo, rho, delta, minIdx ] = densityCluster( X, ...
                                                    d_c, varargin )
%
% DENSITY_CLUSTER -- clustering by fast search and find of the density
%                    peaks
%
% SYNTAX
%
%   [ LABEL, IDXCENTER, BGMASK, DENSITY, DISTANCE ] = ...
%       DENSITYCLUSTER( X, D_C, VARARGIN )
%       
%   [ ... ] = DENSITYCLUSTER( ..., 'nClusters', NCLUSTERS )
%   [ ... ] = DENSITYCLUSTER( ..., 'cutoff', CUTOFF )
%
% INPUT
%
%   X           Feature vectors                         [N-by-D]
%   D_C         Cut-off distance                        [scalar]
%
% OPTIONAL
% 
%   'nClusters' Number of clusters (empty to search)    [scalar]
%               {default: []}
%   'cutoff'    Cut-off type used                       [string]
%                * 'hard': hard cutoff
%                * 'soft': gaussian weighted
%               {default: 'soft'}
% 
% OUTPUT
%
%   LABEL       Cluster labels                          [1-by-N]
%   IDXCENTER   Cluster center indices                  [1-by-C]
%   HALO        Border point mask                       [1-by-N; logical]
%   DENSITY     local density estimate                  [1-by-N]
%   DISTANCE    distance to higher density              [1-by-N]
%
% DESCRIPTION
%
%   [ LABEL, IDXCENTER, BGMASK ] = DENSITYCLUSTER ( X, D_C, VARARGIN ) 
%   returns the result via density-based clustering by fast search and find
%   of density peaks.
%
% REFERENCES
%   
%   [1] Rodriguez, A. and Laio, A., 2014. Clustering by fast search
%       and find of density peaks. Science, 344(6191), pp.1492-1496.
% 
%   [2] Mehmood, Rashid, et al. "Clustering by fast search and merge
%       of local density peaks for gene expression microarray data."
%       Scientific reports 7 (2017): 45602.
%
%   http://people.sissa.it/~laio/Research/Res_clustering.php
%

    %% OPTIONAL PARAMETERS
        
    opt.nClusters = [];
    opt.cutoff    = 'soft';
    
    % parsing the optional arguments
    opt = parseOptArgs( opt, varargin{:} );
    

    %% FIND DENSITY AND DISTANCE TO HIGHER DENSITY
    if(issparse(X))
    % pairwise distances
        D = X;
    else
        D = pdist2( X, X, 'euclidean');
    end
    % get decision graph
    [rho, delta, minIdx, perm] = genDecisionGraph( D, d_c, opt.cutoff );
    
    
    %% LOCATE THE CLUSTER CENTERS BY THRESHOLDING
    
    if isempty( opt.nClusters )
    
      % find cluster centers on the Decision Graph [2]
      idxCenter =  find( and( delta > d_c, ...
                              rho > mean( rho ) ) );
    
    elseif opt.nClusters == 0
      
      % find cluster centers on the Decision Graph [2]
      idxCenter =  find( delta > d_c );
      
    else
      
      % find gamma values
      gamma = rho .* delta;
      
      % get top [opt.nClusters] values
      [~,ii]    = sort( gamma, 'descend' );
      idxCenter = ii( 1:opt.nClusters );
      
    end
      
    
    %% FIND THE MEMBERSHIP OF OTHER POINTS

    % assing each point to closest cluster
    clusterID = findClusterID( minIdx, idxCenter, perm );
    
    
    %% FIND THE POINTS WITH UNCERTAINTY
    
    % get halo mask
    halo = findHalo( D, rho, clusterID, d_c );
    
    
end


%%------------------------------------------------------------
%
% AUTHORS
%
%   Tiancheng Liu                       tcliu@cs.duke.edu
%   Dimitris Floros                      fcdimitr@auth.gr
%
% VERSION
%
%   0.3 (June 07, 2018) - Dimitris
%
% CHANGELOG
%
%
% VERSION
%   0.4 - March 15, 2023
%   
%
% CHANGELOG
%   0.4 (March 15, 2023) - Xenofon
%       *added sparse functionality 
%
%   0.3 (Jun 07, 2018) - Dimitris
%       * added extra optional command to select hard or soft
%         cutoff
% 
%   0.2 (Apr 27, 2018) - Dimitris
%       * broken into sub-functions
%       * fixed minor issues
%       * removed optional inputs
% 
%   0.1 (Oct 12, 2017) - Tiancheng
%       * initial implementation
%
% ------------------------------------------------------------
