%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to compute the likelihood
% p(t2 | t1, b_j1, b_j2)
% by Anirudh Vemula, Jul 20, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [likelihood l_vector] = DP_traj_likelihood_dep(sparseGP_x_2, ...
                                                  sparseGP_y_2, ...
                                                  traj_2, traj_1)


% Lengths of the trajectories
n1 = length(traj_1.x);
n2 = length(traj_2.x);

% Calculate the conditional likelihood 
[mu_x_vec_2, var_x_matrix_2] = sparseGP_predict(sparseGP_x_2, traj_2.x, ...
                                                traj_2.y);
[mu_y_vec_2, var_y_matrix_2] = sparseGP_predict(sparseGP_y_2, traj_2.x, ...
                                                traj_2.y);

var_x_vec_2 = diag(var_x_matrix_2);
var_y_vec_2 = diag(var_y_matrix_2);

% PARAMETERS
h = 0.5;
alpha = 1;

likelihood = 0;
currentTime = 0;
for i=1:n2
   
    currentTime = currentTime + traj_2.dt(i);
    
    mu_x = mu_x_vec_2(i);
    var_x = var_x_vec_2(i);
    mu_y = mu_y_vec_2(i);
    var_y = var_y_vec_2(i);
    
    x_diff = traj_2.dx_dt(i) - mu_x;
    y_diff = traj_2.dy_dt(i) - mu_y;
    
    lx = -0.5 * log(var_x) - 0.5 * x_diff^2 / var_x;
    ly = -0.5 * log(var_y) - 0.5 * y_diff^2 / var_y;
    
    likelihood = likelihood + lx + ly;
    
    if (i < n2)
        % Now for the conditional term
        % Prediction for the next step
        x_new = traj_2.x(i) + mu_x * traj_2.dt(i);
        y_new = traj_2.y(i) + mu_y * traj_2.dt(i);
        
        % Trajectory 1 position
        [x_1, y_1] = findTrajPosition(traj_1, currentTime);
        
        % weight it with how close it is to the traj_1
        % The term we use is (1 - alpha * exp((-1/2h^2) * |t2(x,y) - t1(x,y)|))
        lc = log (1 - alpha * exp( (-1/(2*h^2) * norm([x_new y_new] - ...
                                                      [x_1 y_1]) ) ) );        
        likelihood = likelihood + lc;
    end
    
    %    [lx, ly, lc]
end

% Normalization on the basis of traj length
likelihood = likelihood;% / length(traj_2.x);

end