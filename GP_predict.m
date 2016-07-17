%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 16.322 Stochastic Estimation
% Gaussian Process
% prediction using Gaussian kernal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [mu, var_f] = GP_predict(x_obs, y_obs, speed, x_query, y_query, hyperparam)

% compute C
C = Gaussian_kernel(x_obs, y_obs, x_obs, y_obs, hyperparam);
%C
%C_inv = inv(C);
%x_speed
%y_speed
alpha = C\speed;
mu = zeros(length(x_query),1);
var_f = zeros(length(x_query),length(x_query));

    % plot kernel
    % figure
    % [X,Y] = meshgrid(x_obs,x_obs);
    % contour(X,Y,C, 'ShowText','on')

% compute mean and variance for each data point

% for i = 1: length(x_query)
%     Cxx = Gaussian_kernel(x_query(i), y_query(i), x_query(i), y_query(i), hyperparam,0);
%     kx = Gaussian_kernel(x_obs, y_obs, x_query(i), y_query(i), hyperparam,0);
%     % GP update
%     % consider bias
%     mu(i) = kx'* alpha ;
%     var_f(i) = Cxx - kx' * (C\kx);
% end
Kxx = Gaussian_kernel(x_query, y_query, x_query, y_query, hyperparam,1);
Kxx_prime = Gaussian_kernel(x_query, y_query, x_obs, y_obs, hyperparam,1);
mu = Kxx_prime * alpha;
var_f = Kxx - Kxx_prime * inv(C) * Kxx_prime';

%C
%C_inv
% n = length(x_obs);
% log_like = -n/2*log(2*pi)-1/2*log(det(C))-1/2*y_obs*C_inv*y_obs';

end