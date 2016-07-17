function plotTrajs(trajs, string_in)

n_trajs = trajs.n_traj;
colors = {'r','b','g','y','m','c','k'};

figure
hold on
num_group_per_subplot = 7;
num_cluster = max(trajs.cluster(:,end));
% trajs with cluster assignment 0 are singleton clusters
ifRemove = 1-min(trajs.cluster(:,end)); 
num_figure = ceil(num_cluster / num_group_per_subplot ) + ifRemove;
row_num = ceil(num_figure / 2);
col_num = ceil(num_figure /row_num);

%figure
for k = 1:num_figure
    subplot(row_num,col_num, k);
    low_bnd = (k-1) * num_group_per_subplot ;
    high_bnd = k * num_group_per_subplot ;
    for i=1:n_trajs
        clus_assign = trajs.cluster(i, end);
        if (clus_assign > low_bnd) && (clus_assign <= high_bnd)
            c = colors{mod(trajs.cluster(i,end),7)+1};
            plot(trajs.data(i).x, trajs.data(i).y, c);
            hold on
            plot(trajs.data(i).x(1), trajs.data(i).y(1), strcat(c,'*'));
        end
    end
    
    if k == 1
        if (nargin == 2)
            title(sprintf('%s - found %d clusters, total %d trajs', string_in, num_cluster, n_trajs));
        else 
            title(sprintf('found %d clusters, total %d trajs', num_cluster, n_trajs));
        end
    end
    num_cluster = max(trajs.cluster(:,end));
    axis equal
end

if ifRemove == 1
   subplot(row_num,col_num,num_figure); 
   for i=1:n_trajs
        clus_assign = trajs.cluster(i, end);
        if (clus_assign == 0) 
            c = 'k';
            plot(trajs.data(i).x, trajs.data(i).y, c);
            hold on
            plot(trajs.data(i).x(1), trajs.data(i).y(1), 'k*');
            
        end
   end
   title('singleton clusters, not included in GPs');
end


    set(gcf, 'Position', [100, 100, 300*row_num, 300*col_num]);
end

function sat_value = saturate(l_min, l_max, value)
    sat_value = min(l_max, max(l_min, value));
end