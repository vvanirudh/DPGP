function trajs = generateTrajsWithAvoidanceAndTau(n_traj, n_points, ...
                                                  pLimit, speed, ...
                                                  sigma_noise)

    func = {@linear_func};
    
    num_func = length(func);
    
    traj_counter = 1;
    
    threshold = 0.5;
    
    trajs = struct;
    trajs.data = struct;
    flag = 1;
    
    while (traj_counter <= n_traj*2)       
        speed_in = randn(1) * (0.00 * speed) + speed;
        
        [x1, y1, x2, y2, dt1, dt2] = generateTwoAgentTraj(n_points, pLimit, ...
                                                          speed, ...
                                                          threshold, ...
                                                          func{1});
        
        [t1, t2] = processTrajWithTau(x1, y1, dt1, x2, y2, dt2, ...
                                      sigma_noise);
        
        trajs.data(traj_counter).x = t1.x;
        trajs.data(traj_counter).y = t1.y;
        trajs.data(traj_counter).dx_dt = t1.dx_dt;
        trajs.data(traj_counter).dy_dt = t1.dy_dt;
        trajs.data(traj_counter).dt = dt1;
        trajs.data(traj_counter).speed = speed_in;
        trajs.data(traj_counter).tau = t1.tau;
        
        trajs.data(traj_counter+1).x = t2.x;
        trajs.data(traj_counter+1).y = t2.y;
        trajs.data(traj_counter+1).dx_dt = t2.dx_dt;
        trajs.data(traj_counter+1).dy_dt = t2.dy_dt;
        trajs.data(traj_counter+1).dt = dt2;
        trajs.data(traj_counter+1).speed = speed_in;
        trajs.data(traj_counter+1).tau = t2.tau;
        
        traj_counter = traj_counter + 2;
        flag = 1 - flag;        
    end
    
    trajs.n_traj = n_traj*2;
    trajs.cluster = ones(n_traj*2,1);
    trajs.sweep_count = 1;
    trajs.n_clus = 1;
    
end
