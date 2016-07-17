function likelihood = DP_traj_likelihood(sparseGP_x, sparseGP_y, traj)

[mu_x, var_x] = sparseGP_predict(sparseGP_x, traj.x, traj.y);
[mu_y, var_y] = sparseGP_predict(sparseGP_y, traj.x, traj.y);

x_diff = traj.dx_dt - mu_x;
y_diff = traj.dy_dt - mu_y;

 lx = -0.5 * log(abs(det(var_x))) - 0.5 * x_diff' * (inv(var_x) * x_diff);
% mu_x
% var_x
% -0.5 * log(abs(det(var_x)))
% - 0.5 * x_diff' * (var_x * x_diff)
 ly = -0.5 * log(abs(det(var_y))) - 0.5 * y_diff' * (inv(var_y) * y_diff);
%lx = - 0.5 * x_diff' * (inv(var_x) * x_diff);
%ly = - 0.5 * y_diff' * (inv(var_y) * y_diff);
likelihood = lx + ly;

% normalization based on traj length
likelihood = likelihood/(length(traj.x));
end