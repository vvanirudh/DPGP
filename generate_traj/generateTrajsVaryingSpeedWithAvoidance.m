function trajs = generateTrajsVaryingSpeedWithAvoidance(n_traj, n_points, pLimit, speed, sigma_noise)
% all behavior patterns
%func = {@sinusoidal, @quadratic, @cubic, @linear_func};
%func = {@linear_func, @linear_func2};
    func = {@linear_func};
    %    func = {@linear_func, @quadratic};
    num_func = length(func);
    traj_counter = 1;
    
    threshold = 0.5;
    
    % initialize trajs structure
    % each column corresponds to the trajectory of one agent
    trajs = struct;
    
    trajs.data = struct;
    flag = 1;
    while (traj_counter <= n_traj*2)
        speed_in = randn(1) * (0.00 * speed) + speed;
        
        [x1, y1, x2, y2, dt1, dt2] = generateTwoAgentTraj(n_points, pLimit, ...
                                                          speed, ...
                                                          threshold, ...
                                                          func{1});
        
        
        [x1, y1, dx1_dt, dy1_dt] = processTraj(x1, y1, dt1, ...
                                               sigma_noise);
        [x2, y2, dx2_dt, dy2_dt] = processTraj(x2, y2, dt2, ...
                                               sigma_noise);
        
        
        trajs.data(traj_counter).x = x1;
        trajs.data(traj_counter).y = y1;
        trajs.data(traj_counter).dx_dt = dx1_dt;
        trajs.data(traj_counter).dy_dt = dy1_dt;
        trajs.data(traj_counter).dt = dt1;
        trajs.data(traj_counter).speed = speed_in;
        
        trajs.data(traj_counter+1).x = x2;
        trajs.data(traj_counter+1).y = y2;
        trajs.data(traj_counter+1).dx_dt = dx2_dt;
        trajs.data(traj_counter+1).dy_dt = dy2_dt;
        trajs.data(traj_counter+1).dt = dt2;
        trajs.data(traj_counter+1).speed = speed_in;                
        
        traj_counter = traj_counter + 2;
        flag = 1 - flag;
    end
    
    % initilization
    % every trajectory belongs to the same cluster
    trajs.n_traj = n_traj*2;
    trajs.cluster = ones(n_traj*2,1);
    trajs.sweep_count = 1;
    trajs.n_clus = 1;
    %debugging
    %plotTraj( trajs.data(1).x, trajs.data(1).y, trajs.data(1).dx_dt, trajs.data(1).dy_dt );
    %plotTraj( trajs.data(2).x, trajs.data(2).y, trajs.data(2).dx_dt, trajs.data(2).dy_dt );
    %plotTraj( trajs.data(3).x, trajs.data(3).y, trajs.data(3).dx_dt, trajs.data(3).dy_dt );
    %plotTraj( trajs.data(4).x, trajs.data(4).y, trajs.data(4).dx_dt, trajs.data(4).dy_dt );
end