%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to predict trajectories based
% on best fitting cluster
% by Anirudh Vemula, Jul 28, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function trajs_c = predictTrajs(trajs, num_steps)

global sparseGPs

n_traj = trajs.n_traj;

traj_counter = 1;

trajs_c = struct;
trajs_c.data = struct;
trajs_c.n_traj = n_traj;

while traj_counter <= n_traj
   
    traj = trajs.data(traj_counter);
    clus = trajs.cluster(traj_counter, end);
    
    sparseGP_x = sparseGPs(clus).sparseGP_x;
    sparseGP_y = sparseGPs(clus).sparseGP_y;
    
    pos = [traj.x(end), traj.y(end)];
    pos_var_in = [0.5, 0; 0, 0.5];
    
    tmp_dt = ones(num_steps, 1);
    speed = traj.speed;
    
    pred_traj = sparseGP_trajPredict(sparseGP_x, sparseGP_y, pos, ...
                                     pos_var_in, tmp_dt, num_steps, ...
                                     speed);
    
    predPlot(pred_traj, sparseGP_x, sparseGP_y);
    
    traj_counter = traj_counter + 1;
end

end