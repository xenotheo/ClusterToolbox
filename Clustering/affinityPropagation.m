function [labels, exemplars] = affinityPropagation(X, varargin)
% 
% Affinity Propagation matlab implementation
%         The function takes the N dimmensional points as inputs and
%         returns the Exemplars as well as the label for each point
%   
% SYNTAX
%
%   [ L, E ] = affinityPropagation( X )
%   [ ...  ] = affinityPropagation( X, 'damping', 0.5 )
%
% INPUT
% 
%   X           The dataset of N dimmensional points            [N-by-d]
%   
% OPTIONAL
% 
%   DAMPING     The damping factor for numerical stabilization  [scalar]
%               Authors advised to choose a damping factor 
%               within the range of 0.5 to 1.
%               {default: 0.5}
%   PREFERENCE  The diagonal values of Similarity Matrix        [N-by-1], 
%               {default: 0}                                    [scalar]

%   MAX_ITER    Maximum number of iterations                   [Integer] 
%               {default: 200}
%
%   CONV_ITER   Affinity Propagation will be remarked as       [Integer]
%               converged if within CONV_ITER number of iteratio 
%               the cluster boundaries remain unchanged.
%               {default: 15}
% 
% OUTPUT
% 
%   LABELS      An array contaning the labels for each point    [N-vector]
%               The structure contains the following members:   [N-vector]
%   
%   EXEMPLARS   An array contaning the exemplars of the clusters[Clust-by-1] 
% 
% EXAMPLE
%   load AffinityPoints.mat
%   
%   [L,E] = affinityPropagation(X,'preference', - 50);
%   nC = length(unique(L));
%   figure
%   scatter(X(:,1) , X(:,2) , [] , L);
%   colormap(hsv(nC))
%   hold on
%   scatter(X(E,1) , X(E,2) , [] , L(E),'filled','k')
%   title('Affinity Propagation Clustering Result')
%   hold off

%% PARSE OPTIONAL INPUT
opt = parseOptArgs( varargin{:} );
    
n = size(X,1);

%% CALCULATE SIMILARITY
S = -(pdist2(X,X).^2); 
S(1:n+1:end) = opt.preference;
A = zeros(n,n);
R = zeros(n,n);
eC = zeros(n ,opt.conv_iter);

%% LOOP UNTIL CONVERGENCE OR MAX ITER REACHED
for itter = 1:opt.max_iter
    R = updateR(opt.damping,A,S,R,n); 
    A = updateA(opt.damping,A,R,n);
    exemplars = (diag(R)+diag(A))>0;
    K = sum(exemplars);
    eC(:,mod(itter, opt.conv_iter)+1) = exemplars; 
    if (itter >= opt.conv_iter)
        se  = sum(eC , 2);
        uncov = ( sum((se == opt.conv_iter) + (se ==0))  ~= n);
        if (~(uncov) && (K>0))
           break
        end
    end
end

if K>0
    [labels,~] = knnsearch(X(exemplars,:),X);
    exemplars = find(exemplars);
else
    warning('Affinity Propagation did not converge. There are no exemplars');
    warning('You might want to lower damping or setting a higher Value for Preference');
    labels = -1*ones(n,1);
end
    
        

    

%% LOCAL FUNCTION: RESPONSIBILITY MATRIX UPDATE
function [Rl] = updateR(damping,Al,Sl,Rl,n)
        V = Al + Sl;
    
        %First Max
        [D1 , I1] = max(V,[],2);
        
        %Replace the maximum values with -inf
        V(sub2ind([n,n],(1:n)',I1)) = -inf;
        
        %Second Max
        D2 = max(V,[],2);
        V = Sl - D1;
        V(sub2ind([n,n],(1:n)',I1)) = Sl(sub2ind([n,n],(1:n)',I1)) - D2;
        V = V * (1 - damping);

        Rl = Rl * damping + V;
end

%% LOCAL FUNCTION: AVAILABILITY MATRIX UPDATE
function [Al] = updateA(damping,Al,Rl,n)
        a = max(Rl,0);
        a(1:n+1:end) = Rl(1:n+1:end);
        a = a - sum(a,1);
        dA = diag(a);
        a(a<0) = 0;
        a(1:n+1:end) = dA;        
        Al = (Al * damping) - (1-damping)*a;
        
end

%% LOCAL FUNCTION: OPTIONAL ARGUMENT PARSING
%
function opt = parseOptArgs(varargin)

% default values
dflt.max_iter      = 200;
dflt.damping       = 0.5;
dflt.conv_iter     = 15;
dflt.distPreComp   = false;
dflt.preference    = 0;
dflt.metric        = 'euclidean';
dflt.verbose       = false; 
% inputParser initialization
ip = inputParser;

ip.CaseSensitive   = false;
ip.KeepUnmatched   = false;
ip.PartialMatching = true;
ip.StructExpand    = true;

% set parameters & default values
argNames = fieldnames( dflt );
for i = 1 : length(argNames)
  addParameter( ip, argNames{i}, dflt.(argNames{i}) );
end

% parse input
parse( ip, varargin{:} );
opt = ip.Results;

% use default values for any options set to []
for i = 1 : length(argNames)
  if isempty( opt.(argNames{i}) )
    opt.(argNames{i}) = dflt.(argNames{i});
  end
end
end


end
