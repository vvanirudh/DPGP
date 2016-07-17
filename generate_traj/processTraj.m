function [x, y, dx_dt, dy_dt] = processTraj(x, y, dt, sigma_noise)

    n = length(x);


    
    dx_dt = zeros(n,1);
    dy_dt = zeros(n,1);


    for i=1:n-1
        dx_dt(i) = (x(i+1) - x(i)) / dt(i); 
        dy_dt(i) = (y(i+1) - y(i)) / dt(i);
    end

    
    dx_dt(end,1) = dx_dt(end-1,1);
    dy_dt(end,1) = dy_dt(end-1,1);
    
        
    %rng(1);
    dx_noise = sigma_noise * randn(n,1);
    dy_noise = sigma_noise * randn(n,1);
    
    dx_dt = dx_dt + dx_noise;
    dy_dt = dy_dt + dy_noise;
    
    x(2:n) = x(1:n-1) + dx_dt(1:n-1) .* dt(1:n-1);
    y(2:n) = y(1:n-1) + dy_dt(1:n-1) .* dt(1:n-1);

end