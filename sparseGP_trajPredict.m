function pred_traj = sparseGP_trajPredict(sparseGP_x, sparseGP_y, pos, pos_var_in, dt, num_steps, speed)

if length(dt) ~= num_steps
    dt = dt * ones(1,num_steps);
end

pos_var = pos_var_in;
pred_traj = zeros(num_steps+1,4);
pred_traj(1, 1:2) = pos;
pred_traj(1, 3:4) = [pos_var(1,1), pos_var(2,2)];

for i = 1:num_steps
   [mu_vx, var_vx] = sparseGP_predict_distri(sparseGP_x, pred_traj(i, 1:2), pos_var);
   [mu_vy, var_vy] = sparseGP_predict_distri(sparseGP_y, pred_traj(i, 1:2), pos_var);
   pred_traj(i+1,1:2) = pred_traj(i,1:2) + dt(i) * [mu_vx, mu_vy] * speed/sqrt(mu_vx^2 + mu_vy^2);
   pred_traj(i+1,3:4) = pred_traj(i,3:4) + dt(i)^2 * [var_vx, var_vy];
   pos_var = diag(pred_traj(i+1,3:4));
end

%fprintf(sprintf('speed: %f',speed))
%dt

end