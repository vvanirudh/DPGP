function pred_traj = linearTrajPred(pos, pos_var_in, dt, num_steps, speed_x, speed_y)

if length(dt) ~= num_steps
    dt = dt * ones(1,num_steps);
end

pos_var = pos_var_in;
pred_traj = zeros(num_steps+1,4);
pred_traj(1, 1:2) = pos;
pred_traj(1, 3:4) = [pos_var(1,1), pos_var(2,2)];

var_vx = 0.2; 
var_vy = 0.2;

for i = 1:num_steps
   pred_traj(i+1,1:2) = pred_traj(i,1:2) + dt(i) * [speed_x, speed_y];
   pred_traj(i+1,3:4) = pred_traj(i,3:4) + dt(i)^2 * [var_vx, var_vy];
end

end
