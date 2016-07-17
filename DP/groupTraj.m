function count = groupTraj(sweep_count)
%fprintf('in groupTraj function, sweep number: %d', sweep_count);
global trajs sparseGPs
% preprocessing: remove empty clusters
[c, ia, ic] = unique(trajs.cluster(:,sweep_count));
[c_sort, c_ind] = sort(c);
trajs.cluster(:,sweep_count) = c_ind(ic);
trajs.cluster(:,sweep_count) = reorder_vec(trajs.cluster(:,sweep_count)')';
trajs.n_clus = max(trajs.cluster(:,sweep_count));

count = ones(trajs.n_clus,1);
% empty sparseGPs structure
sparseGPs = {};
    for i = 1:trajs.n_clus
        indices = find(trajs.cluster(:,sweep_count) == i);
        count(i) = length(indices);
        % concacenate all trajectories belonging to the same cluster
        % into one GP
        sparseGPs(i).data.x = [];
        sparseGPs(i).data.y = [];
        sparseGPs(i).data.dx_dt = [];
        sparseGPs(i).data.dy_dt = [];
        sparseGPs(i).count = count(i);
        for j = 1:length(indices)
            sparseGPs(i).data.x = [sparseGPs(i).data.x; trajs.data(indices(j)).x];
            sparseGPs(i).data.y = [sparseGPs(i).data.y; trajs.data(indices(j)).y];
            sparseGPs(i).data.dx_dt = [sparseGPs(i).data.dx_dt; trajs.data(indices(j)).dx_dt];
            sparseGPs(i).data.dy_dt = [sparseGPs(i).data.dy_dt; trajs.data(indices(j)).dy_dt];
            sparseGPs(i).trajs(j).dx_dt = trajs.data(indices(j)).dx_dt;
            sparseGPs(i).trajs(j).dy_dt = trajs.data(indices(j)).dy_dt;
        end        
    end


end

function vec_p = reorder_vec(vec)
    n = length(vec);
    vec_p = 9999 * ones(1,n);
    num = 1;
    for i = 1:n
       if vec_p(i) == 9999;
          ind = find(vec == vec(i));
          vec_p(ind) = num;
          num = num + 1;
          
       end
    end   
end