function l = self_likelihood_indep(traj, hyperparam)

% fit a GP
sparseGP_x = build_sparseGP(traj.x, traj.y, traj.dx_dt, hyperparam,0);
sparseGP_y = build_sparseGP(traj.x, traj.y, traj.dy_dt, hyperparam,0);

% compute likelihood
l = DP_traj_likelihood_indep(sparseGP_x, sparseGP_y, traj);

end