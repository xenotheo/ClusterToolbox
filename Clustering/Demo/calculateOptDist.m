function [optDist] = calculateOptDist(X , minPts)

[~,X] = knnsearch(X,X,'k',minPts);
X = sort(X(:,minPts));
SD = smooth(X, max(ceil(size(X,1)/34) , 13), 'lowess');
slide = arrayfun(@(x,y,z) x - y , SD(2:end-1) , SD(1:end-2));
curv  = arrayfun(@(x,y,z) x +y -2*z,SD(3:end) , SD(1:end-2) , SD(2:end-1));
idx = (2 : size(X,1)-1)';
taylor1 = SD(idx) + (idx' - idx).*slide(idx - 1);
taylor2 = SD(idx) + (idx' - idx).*slide(idx - 1) + 0.5*((idx' - idx).^2).*curv(idx - 1);

error = sum(0.5*(abs(taylor1 - SD(idx))) + 0.5*abs(taylor2 - SD(idx)) ,1);
[~,I] = min(error);
optDist = X(I+1);
end