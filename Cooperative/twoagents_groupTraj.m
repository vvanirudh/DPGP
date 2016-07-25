%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to group two agent trajectories
% by Anirudh Vemula, Jul 19, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function count = twoagents_groupTraj(sweep_count)

global trajs sparseGPs

% preprocessing: remove empty clusters
[c, ia, ic] = unique(trajs.cluster(:, :, sweep_count));
[c_sort, c_ind] = sort(c);
trajs.cluster(:, :, sweep_count) = reshape(c_ind(ic), ...
                                           size(trajs.cluster,1), ...
                                           size(trajs.cluster, 2));

trajs.cluster(:, :, sweep_count) = reorder_matrix(trajs.cluster(:, ...
                                                  :, sweep_count));
trajs.n_clus = max(max(trajs.cluster(:,:,sweep_count)));

count = ones(trajs.n_clus, 1);
% Empty sparseGPs structure
sparseGPs = {};

for i=1:trajs.n_clus
   
    [indices_row indices_column] = find(trajs.cluster(:,:,sweep_count) == i);
    count(i) = length(indices_row);
    % concatenate all trajectories belonging to the same cluster
    % into one GP
    sparseGPs(i).data.x = [];
    sparseGPs(i).data.y = [];
    sparseGPs(i).data.dx_dt = [];
    sparseGPs(i).data.dy_dt = [];
    sparseGPs(i).count = count(i);
    
    for j=1:length(indices_row)
       if indices_column(j)==1 
           % Trajectory of 1st agent
           sparseGPs(i).data.x =  [sparseGPs(i).data.x; ...
                               trajs.data1(indices_row(j)).x];
           sparseGPs(i).data.y =  [sparseGPs(i).data.y; ...
                               trajs.data1(indices_row(j)).y];
           sparseGPs(i).data.dx_dt =  [sparseGPs(i).data.dx_dt; ...
                               trajs.data1(indices_row(j)).dx_dt];
           sparseGPs(i).data.dy_dt =  [sparseGPs(i).data.dy_dt; ...
                               trajs.data1(indices_row(j)).dy_dt];
           sparseGPs(i).trajs(j).dx_dt = ...
               trajs.data1(indices_row(j)).dx_dt;
           sparseGPs(i).trajs(j).dy_dt = ...
               trajs.data1(indices_row(j)).dy_dt;
       else
           % Trajectory of 2nd agent
           sparseGPs(i).data.x =  [sparseGPs(i).data.x; ...
                               trajs.data2(indices_row(j)).x];
           sparseGPs(i).data.y =  [sparseGPs(i).data.y; ...
                               trajs.data2(indices_row(j)).y];
           sparseGPs(i).data.dx_dt =  [sparseGPs(i).data.dx_dt; ...
                               trajs.data2(indices_row(j)).dx_dt];
           sparseGPs(i).data.dy_dt =  [sparseGPs(i).data.dy_dt; ...
                               trajs.data2(indices_row(j)).dy_dt];
           sparseGPs(i).trajs(j).dx_dt = ...
               trajs.data2(indices_row(j)).dx_dt;
           sparseGPs(i).trajs(j).dy_dt = ...
               trajs.data2(indices_row(j)).dy_dt;
       end
    end    
end

end

function mat_p = reorder_matrix(mat)   
    [m, n] = size(mat);
    mat_p = 9999 * ones(m, n);
    num = 1;
    for i=1:m
        for j=1:n
            if mat_p(i,j) == 9999
                ind = find(mat == mat(i,j));
                mat_p(ind) = num;
                num = num + 1;
            end
        end
    end    
end