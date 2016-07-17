function [mu, var]= sparseGP_predict(sparseGP, x_query, y_query)

% compute C
mu = zeros(length(x_query),1);
hyperparam = sparseGP.hyperparam;

    % plot kernel
    % figure
    % [X,Y] = meshgrid(x_obs,x_obs);
    % contour(X,Y,C, 'ShowText','on')

% compute mean and variance for each data point
kx = Gaussian_kernel(x_query, y_query, sparseGP.BV_x, sparseGP.BV_y, hyperparam,0);
kxx = Gaussian_kernel(x_query, y_query, x_query, y_query, hyperparam,1);
mu = kx * sparseGP.alpha;
var = abs(kxx + kx * sparseGP.C * kx');

% to account for variance due to shift
x_trans = [-0.5, 0, 0.5];
y_trans = [-0.5, 0, 0.5];

x_query_aug = [x_query+x_trans(1); x_query+x_trans(2); x_query+x_trans(3); ...
            x_query+x_trans(1); x_query+x_trans(3); ...
            x_query+x_trans(1); x_query+x_trans(2); x_query+x_trans(3)];
        
y_query_aug = [y_query+y_trans(1); y_query+y_trans(1); y_query+y_trans(1); ...
            y_query+y_trans(2); y_query+y_trans(2); ...
            y_query+y_trans(3); y_query+x_trans(3); y_query+y_trans(3)];
        
kx_aug = Gaussian_kernel(x_query_aug, y_query_aug, sparseGP.BV_x, sparseGP.BV_y, hyperparam,0);
mu_aug = kx_aug * sparseGP.alpha;
% modify diagonal of variance matrix
t = length(x_query);

%for i = 1:t
%  var(1 + (i-1)*(t+1)) = min(var(1 + (i-1)*(t+1)),hyperparam(4)^2) +...
%      max(abs(mu_aug(i:t:end) - mu(i)))^2;
%end

% % for trajectory translation
% if (size(x_query) == [1,1])
%     x_trans = [-0.5, 0, 0.5];
%     y_trans = [-0.5, 0, 0.5];
%     mu_around = mu * ones(3,3);
%     for i = 1:3
%         for j = 1:3
%             if (i == 2 && j == 2)
%                 continue;
%             else
%                 tmp = Gaussian_kernel(x_query+x_trans(i), y_query+y_trans(j), sparseGP.BV_x, sparseGP.BV_y, hyperparam,0);
%                 mu_around(i,j) = tmp * sparseGP.alpha;
%             end
%         end
%     end
%     var = min(var, 2 * hyperparam(4)^2) + max(abs(reshape(mu_around,[9,1]) - mu))^2;
%     var =  hyperparam(4)^2 + max(abs(reshape(mu_around,[9,1]) - mu))^2;
% end

% if var < 0
%    fprintf(sprintf('var = %f, must be an error here!\n',var)); 
%    kx * sparseGP.C * kx'
%    kxx
% end
%for i = 1: length(x_query)
%    kx = Gaussian_kernel(sparseGP.BV_x, sparseGP.BV_y, x_query(i), y_query(i), hyperparam,0);
%    % GP update
%    mu(i) = kx * sparseGP.alpha;
%end


%C
%C_inv
% n = length(x_obs);
% log_like = -n/2*log(2*pi)-1/2*log(det(C))-1/2*y_obs*C_inv*y_obs';

end