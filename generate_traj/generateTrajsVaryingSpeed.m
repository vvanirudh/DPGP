function trajs = generateTrajsVaryingSpeed(n_traj, n_points, pLimit, speed, sigma_noise)
    % all behavior patterns
    func = {@sinusoidal, @quadratic, @cubic, @linear_func};
    num_func = length(func);
    traj_counter = 1;
    % initialize trajs structure
    % each column corresponds to the trajectory of one agent
    trajs = struct;
    
    trajs.data = struct;
    while (traj_counter <= n_traj)
        speed_in = randn(1) * (0.00 * speed) + speed;
        [x, y, dt] = generateTraj(n_points, pLimit, speed_in, func{mod(traj_counter,num_func)+1});
        [x, y, dx_dt, dy_dt] = processTraj(x, y, dt,  sigma_noise);
        
        trajs.data(traj_counter).x = x;
        trajs.data(traj_counter).y = y;
        trajs.data(traj_counter).dx_dt = dx_dt;
        trajs.data(traj_counter).dy_dt = dy_dt;
        trajs.data(traj_counter).dt = dt;
        trajs.data(traj_counter).speed = speed_in;
        traj_counter = traj_counter + 1;
    end
    
    % initilization
    % every trajectory belongs to the same cluster
    trajs.n_traj = n_traj;
    trajs.cluster = ones(n_traj,1);
    trajs.sweep_count = 1;
    trajs.n_clus = 1;
    %debugging
    %plotTraj( trajs.data(1).x, trajs.data(1).y, trajs.data(1).dx_dt, trajs.data(1).dy_dt );
    %plotTraj( trajs.data(2).x, trajs.data(2).y, trajs.data(2).dx_dt, trajs.data(2).dy_dt );
    %plotTraj( trajs.data(3).x, trajs.data(3).y, trajs.data(3).dx_dt, trajs.data(3).dy_dt );
    %plotTraj( trajs.data(4).x, trajs.data(4).y, trajs.data(4).dx_dt, trajs.data(4).dy_dt );
end