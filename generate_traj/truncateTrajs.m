%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to truncate the prediction trajectories 
% so that the model can be tested in prediction tasks
% by Anirudh Vemula, Jul 28, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function trajs_p = truncateTrajs(trajs)

trajs_p = struct;
trajs_p.data = struct;

n_traj = trajs.n_traj;

traj_counter = 1;

while traj_counter <= n_traj
   
    len = length(trajs.data(traj_counter).x);
    
    trunc = floor(len/3);
    
    trajs_p.data(traj_counter).x = trajs.data(traj_counter).x(1: ...
                                                      trunc);
    trajs_p.data(traj_counter).y = trajs.data(traj_counter).y(1: ...
                                                      trunc);
    trajs_p.data(traj_counter).dx_dt = trajs.data(traj_counter).dx_dt(1: ...
                                                      trunc);
    trajs_p.data(traj_counter).dy_dt = trajs.data(traj_counter).dy_dt(1: ...
                                                      trunc);
    trajs_p.data(traj_counter).dt = trajs.data(traj_counter).dt(1: ...
                                                      trunc);
    trajs_p.data(traj_counter).speed = ...
        trajs.data(traj_counter).speed;
    
    traj_counter = traj_counter + 1;
    
end

trajs_p.n_traj = n_traj;

trajs_p.cluster = ones(n_traj, 1);
trajs_p.sweep_count = 1;
trajs_p.n_clus = 1;

end