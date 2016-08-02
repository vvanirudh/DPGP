function [mu, var] = sparseGP_predict_tau(sparseGP, x_query, y_query, ...
                                          tau_query)

    % compute C
    mu = zeros(length(x_query),1);
    hyperparam = sparseGP.hyperparam;

    % plot kernel
    % figure
    % [X,Y] = meshgrid(x_obs,x_obs);
    % contour(X,Y,C, 'ShowText','on')

    % compute mean and variance for each data point
    in_query = [x_query y_query tau_query];
    in_bv = [sparseGP.BV_x sparseGP.BV_y sparseGP.BV_tau];
    kx = Gaussian_kernel_tau(in_query, in_bv, hyperparam,0);
    kxx = Gaussian_kernel_tau(in_query, in_query, hyperparam,1);
    mu = kx * sparseGP.alpha;
    var = abs(kxx + kx * sparseGP.C * kx');

    % to account for variance due to shift
    x_trans = [-0.5, 0, 0.5];
    y_trans = [-0.5, 0, 0.5];
    tau_trans = [-0.5, 0, 0.5];

    x_query_aug = [x_query+x_trans(1); x_query+x_trans(2); x_query+x_trans(3); ...
                   x_query+x_trans(1); x_query+x_trans(3); ...
                   x_query+x_trans(1); x_query+x_trans(2); x_query+x_trans(3)];
    
    y_query_aug = [y_query+y_trans(1); y_query+y_trans(1); y_query+y_trans(1); ...
                   y_query+y_trans(2); y_query+y_trans(2); ...
                   y_query+y_trans(3); y_query+x_trans(3); ...
                   y_query+y_trans(3)];
    
    tau_query_aug = [tau_query+tau_trans(1); tau_query+tau_trans(2); tau_query+tau_trans(3); ...
                   tau_query+tau_trans(1); tau_query+tau_trans(3); ...
                   tau_query+tau_trans(1); tau_query+tau_trans(2); tau_query+tau_trans(3)];
    
    in_query_aug = [x_query_aug y_query_aug tau_query_aug];
    kx_aug = Gaussian_kernel(in_query_aug, in_bv, hyperparam,0);
    mu_aug = kx_aug * sparseGP.alpha;
    % modify diagonal of variance matrix
    t = length(x_query);

end
