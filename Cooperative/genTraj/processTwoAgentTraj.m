%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to process a trajectory
% by Anirudh Vemula, Jul 19, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, y, dx_dt, dy_dt] = processTwoAgentTraj(x, y, dt, ...
                                                  sigma_noise)
    
    % Length of trajectory
    n = length(x);
    
    % Initialize velocity arrays
    dx_dt = zeros(n, 1);
    dy_dt = zeros(n, 1);
    
    % Fill in the velocity arrays
    for i=1:n-1
        dx_dt(i) = (x(i+1) - x(i)) / dt(i);
        dy_dt(i) = (y(i+1) - y(i)) / dt(i);
    end
    
    dx_dt(end, 1) = dx_dt(end-1, 1);
    dy_dt(end, 1) = dy_dt(end-1, 1);
    
    % Noise
    dx_noise = sigma_noise * randn(n, 1);
    dy_noise = sigma_noise * randn(n, 1);
    
    % Add the noise to velocity
    dx_dt = dx_dt + dx_noise;
    dy_dt = dy_dt + dy_noise;
    
    % Process trajectory to include the noise
    x(2:n) = x(1:n-1) + dx_dt(1:n-1) .* dt(1:n-1);
    y(2:n) = y(1:n-1) + dy_dt(1:n-1) .* dt(1:n-1);
    
end