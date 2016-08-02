function likelihood = DP_traj_likelihood_indep_tau(sparseGP_x, ...
                                                  sparseGP_y, traj)

    n = length(traj.x);

    [mu_x_vec, var_x_matrix] = sparseGP_predict_tau(sparseGP_x, traj.x, ...
                                                    traj.y, traj.tau);
    [mu_y_vec, var_y_matrix] = sparseGP_predict_tau(sparseGP_y, traj.x, ...
                                                    traj.y, traj.tau);

    var_x_vec = diag(var_x_matrix);
    var_y_vec = diag(var_y_matrix);

    likelihood = 0;
    %fprintf('in DP_traj_likelihood_indep\n');
    for i = 1:n
        %[mu_x, var_x] = sparseGP_predict(sparseGP_x, traj.x(i), traj.y(i));
        %[mu_y, var_y] = sparseGP_predict(sparseGP_y, traj.x(i), traj.y(i));

        % speed up
        mu_x = mu_x_vec(i);
        var_x = var_x_vec(i);
        mu_y = mu_y_vec(i);
        var_y = var_y_vec(i);

        x_diff = traj.dx_dt(i) - mu_x;
        y_diff = traj.dy_dt(i) - mu_y;
        %fprintf(sprintf('at (%.2f %.2f), predicted:(%.2f %.2f) traj(%.2f %.2f), var_rt(%.5f %.5f)\n',...
        %    traj.x(i), traj.y(i), mu_x, mu_y, traj.dx_dt(i), traj.dy_dt(i),var_x^0.5, var_y^0.5));
        lx = -0.5 * log(var_x) - 0.5 * x_diff^2 / var_x;
        % mu_x
        % var_x
        % -0.5 * log(abs(det(var_x)))
        % - 0.5 * x_diff' * (var_x * x_diff)
        ly = -0.5 * log(var_y) - 0.5 * y_diff^2 /var_y;
        %lx = - 0.5 * x_diff' * (inv(var_x) * x_diff);
        %ly = - 0.5 * y_diff' * (inv(var_y) * y_diff);
        %if (mu_x < 0.05)
        %   lx = 0;
        %end
        %if (mu_y < 0.05)
        %   ly = 0;
        %end
        likelihood = likelihood + lx + ly;
    end
    % normalization based on traj length
    likelihood = likelihood; %/(length(traj.x));
    
    if isnan(likelihood)
        keyboard()
    end
    
end
