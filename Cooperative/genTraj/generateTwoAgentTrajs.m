%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to generate n trajectories of 
% two interacting agents. For now modeling
% cooperative collision avoidance
% by Anirudh Vemula, Jul 18, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function trajs = generateTwoAgentTrajs(n_traj, n_points, pLimit, ...
                                       speed, sigma_noise)

% Functions used to generate trajectories
%func = {@sinusoidal, @quadratic, @cubic, @linear_func};
%func = {@linear_func, @quadratic};
%func = {@linear_func};
func = {@quadratic};

% SET THRESHOLD
threshold = 0.5;

% Number of functions used
num_func = length(func);

traj_counter = 1;

% Initialize the trajs structure
% Each column corresponds to the pair of trajectories of the two agents
trajs = struct;
trajs.data1 = struct;
trajs.data2 = struct;

while (traj_counter <= n_traj)
   
    speed_in = randn(1) * (0.00 * speed) + speed;
    [x1, y1, x2, y2, dt1, dt2] = generateTwoAgentTraj(n_points, ...
                                                      pLimit, speed_in, ...
                                                      threshold, ...
                                                      func{mod(traj_counter, num_func) + 1});
    
    [x1, y1, dx1_dt, dy1_dt] = processTwoAgentTraj(x1, y1, dt1, ...
                                                   sigma_noise);
    [x2, y2, dx2_dt, dy2_dt] = processTwoAgentTraj(x2, y2, dt2, ...
                                                   sigma_noise);
    
    trajs.data1(traj_counter).x = x1;
    trajs.data1(traj_counter).y = y1;
    trajs.data1(traj_counter).dx_dt = dx1_dt;
    trajs.data1(traj_counter).dy_dt = dy1_dt;
    trajs.data1(traj_counter).dt = dt1;
    trajs.data1(traj_counter).speed = speed_in;
    
    trajs.data2(traj_counter).x = x2;
    trajs.data2(traj_counter).y = y2;
    trajs.data2(traj_counter).dx_dt = dx2_dt;
    trajs.data2(traj_counter).dy_dt = dy2_dt;
    trajs.data2(traj_counter).dt = dt2;
    trajs.data2(traj_counter).speed = speed_in;
    
    traj_counter = traj_counter + 1;
end

% initialization
% Every trajectory belongs to the same cluster
trajs.n_traj = n_traj;

% Each pair of traj has 2 cluster assignments
trajs.cluster = ones(n_traj, 2);

trajs.sweep_count = 1;
trajs.n_clus = 1;


% Debugging
%plotTraj( trajs.data1(1).x, trajs.data1(1).y, trajs.data1(1).dx_dt, ...
%          trajs.data1(1).dy_dt);
%plotTraj( trajs.data2(1).x, trajs.data2(1).y, trajs.data2(1).dx_dt, ...
%          trajs.data2(1).dy_dt);

end