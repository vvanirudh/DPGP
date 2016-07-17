function [f, g] = GP_likelihood(length_param)

    global x_obs_opt y_obs_opt speed_opt hyperparam_opt;
    hyperparam_opt(1) = length_param(1);
    hyperparam_opt(2) = length_param(2);
    sigma_input = hyperparam_opt(3);
    sigma_noise = hyperparam_opt(4);
    lx = length_param(1);
    ly = length_param(2);
    
    % compute C
    % hyperparam
    K = Gaussian_kernel(x_obs_opt, y_obs_opt, x_obs_opt, y_obs_opt, hyperparam_opt, 0);
    
    C = K + sigma_noise^2 * eye(length(x_obs_opt));
    C_inv = inv(C); 
    

    % compute likelihood 
    n = length(x_obs_opt);
    %-n/2*log(2*pi);
    %-1/2*log(det(C));
    %-1/2*speed'*C_inv*speed;
    f = -n/2*log(2*pi)-1/2*log(det(C))-1/2*speed_opt'*(C\speed_opt);
    
    % compute gradient
    [X1,X2] = meshgrid(x_obs_opt,x_obs_opt);
    [Y1,Y2] = meshgrid(y_obs_opt,y_obs_opt);
    dk_dlx = C.* ((X1-X2).^2) / lx^3;
    dk_dly = C.* ((Y1-Y2).^2) / ly^3;

    alpha = C \ speed_opt;
    g(1) = 0.5 * trace( (alpha * alpha' - C_inv) * dk_dlx);
    g(2) = 0.5 * trace( (alpha * alpha' - C_inv) * dk_dly);
    f = -f; % minimize the negative of likelihood
    g = -g;
    %[f g l]
end