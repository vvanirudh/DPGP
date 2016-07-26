%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to generate partial trajectories from
% complete trajectories
% by Anirudh Vemula, Jul 26, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function trajs_partial = twoagents_partialTrajs(trajs)

    trajs_partial = struct;
    trajs_partial.data1 = struct;
    trajs_partial.data2 = struct;
    
    n_traj = trajs.n_traj;
    
    traj_counter = 1;
    
    while (traj_counter <= n_traj)
        
        len1 = length(trajs.data1(traj_counter).x);
        len2 = length(trajs.data2(traj_counter).x);
        
        trajs_partial.data1(traj_counter).x = ...
            trajs.data1(traj_counter).x(1:len1/2);
        trajs_partial.data1(traj_counter).y = ...
            trajs.data1(traj_counter).y(1:len1/2);
        trajs_partial.data1(traj_counter).dx_dt = ...
            trajs.data1(traj_counter).dx_dt(1:len1/2);
        trajs_partial.data1(traj_counter).dy_dt = ...
            trajs.data1(traj_counter).dy_dt(1:len1/2);
        trajs_partial.data1(traj_counter).dt = ...
            trajs.data1(traj_counter).dt(1:len1/2);
        trajs_partial.data1(traj_counter).speed = ...
            trajs.data1(traj_counter).speed;
        
        trajs_partial.data2(traj_counter).x = ...
            trajs.data2(traj_counter).x(1:len2/2);
        trajs_partial.data2(traj_counter).y = ...
            trajs.data2(traj_counter).y(1:len2/2);
        trajs_partial.data2(traj_counter).dx_dt = ...
            trajs.data2(traj_counter).dx_dt(1:len2/2);
        trajs_partial.data2(traj_counter).dy_dt = ...
            trajs.data2(traj_counter).dy_dt(1:len2/2);
        trajs_partial.data2(traj_counter).dt = ...
            trajs.data2(traj_counter).dt(1:len2/2);
        trajs_partial.data2(traj_counter).speed = ...
            trajs.data2(traj_counter).speed;
        
        traj_counter = traj_counter + 1;
    end
    
    trajs_partial.n_traj = n_traj;
    
    trajs_partial.cluster = ones(n_traj, 2);
    trajs_partial.sweep_count = 1;
    trajs_partial.n_clus = 1;
    
end
