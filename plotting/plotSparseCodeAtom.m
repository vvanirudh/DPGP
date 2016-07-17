function plotSparseCodeAtom(globalInfo, count, direction_x, direction_y)

num_figure = 1;
if nargin > 2
    num_figure = 2;
    count = direction_x.^2 > 0;
end

max_count = max(max(count));
if max_count > 9
    max_count = max_count / 3;
end

figure
subplot(num_figure,1, 1)
clims = [0, max_count];
imagesc(flipud(count'),clims)
colormap(gray)
axis image
title('occupancy')

x_stride = (globalInfo.bnd(2) - globalInfo.bnd(1)) / 4;
xticklabels = globalInfo.bnd(1):x_stride:globalInfo.bnd(2);
xticks = linspace(1, size(count,1), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);

y_stride = (globalInfo.bnd(4) - globalInfo.bnd(3)) / 4;
yticklabels = flip(globalInfo.bnd(3):y_stride:globalInfo.bnd(4));
yticks = linspace(1, size(count,2), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);

if num_figure == 2
    subplot(num_figure, 1,2);
    [x, y] = meshgrid(globalInfo.x_ticks+0.5 * globalInfo.width, ...
       globalInfo.y_ticks+0.5 * globalInfo.width);
    quiver(x,y,direction_x',direction_y','r');    
    hold on
    % plot grid
    for i = 1:length(globalInfo.x_ticks)
       plot ([globalInfo.x_ticks(i), globalInfo.x_ticks(i)], globalInfo.bnd(3:4), 'k');
    end
    for j = 1:length(globalInfo.y_ticks)
       plot (globalInfo.bnd(1:2), [globalInfo.y_ticks(j), globalInfo.y_ticks(j)], 'k'); 
    end
    
    title('flow direction');
    axis(globalInfo.bnd);
   %axis equal
    set(gcf, 'Position', [1100, 100, 600, 800]); 
end

end