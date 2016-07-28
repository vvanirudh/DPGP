function l = self_likelihood(traj, hyperparam)

% fit a GP
sparseGP_x = build_sparseGP(traj.x, traj.y, traj.dx_dt, hyperparam);
sparseGP_y = build_sparseGP(traj.x, traj.y, traj.dy_dt, hyperparam);

% compute likelihood and divide by 2
l = DP_traj_likelihood(sparseGP_x, sparseGP_y, traj);

end