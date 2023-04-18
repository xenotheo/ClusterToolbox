function []  = reachabilityPlot(reachability , ordering, epsdbscan, X, labels)

% Utility Function to visualize the ordering of points of OPTICS
%

n = size(X,1);
subplot(2,1,1)
hold on
plot(reachability(ordering),'-o','Color','b','MarkerSize',4)

plot(1:n,ones(n,1)*epsdbscan ,'--','Color','k')

title('Reachability Plot')
xlabel('Ordering')
ylabel('Reachability')
legend({'Reachability', ['dbscan Cutoff = ',num2str(epsdbscan)]},'Location','northeast')
hold off
%% Visualize Clusters with value of eps with dbscan cutoff

nC = length(unique(labels));

subplot(2,1,2)
scatter(X(:,1),X(:,2), [], labels)
hold on
scatter(X(labels==-1,1),X(labels==-1,2),[],'k');
colormap(hsv(nC))
title(['Dbscan Cutoff Clustering eps = ',num2str(epsdbscan)])

hold off
end
