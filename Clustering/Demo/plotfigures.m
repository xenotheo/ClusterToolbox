function [] = plotfigures(algorithms,dataset)
    
    for i = 1:size(algorithms,2)
       for j = 1:size(algorithms{i}.results,1)
          figure
          scatter(dataset{j}.Dataset(:,1) , dataset{j}.Dataset(:,2) , [] , algorithms{i}.results{j});
          title([algorithms{i}.title, ' for dataset ' , ...
                 dataset{j}.Description])
          axis image
          colormap(jet(length(unique(algorithms{i}.results{j}))))
       end
        
        
    end

arrangeFigures(36,6,6);
end