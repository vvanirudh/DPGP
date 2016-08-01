function [t1, t2] = processTrajWithTau(x1, y1, dt1, x2, y2, dt2, ...
                                       sigma_noise)
   
    t1 = struct;
    t2 = struct;
    
    n1 = length(x1);
    n2 = length(x2);
    
    dx1_dt = zeros(n, 1);
    dy1_dt = zeros(n, 1);
    tau1 = zeros(n, 1);
    
    dx2_dt = zeros(n, 1);
    dy2_dt = zeros(n, 1);
    tau2 = zeros(n, 1);
    
    for i=1:n1-1
        dx1_dt(i) = (x1(i+1) - x1(i)) / dt1(i);
        dy1_dt(i) = (y1(i+1) - y1(i)) / dt1(i);
    end
    
    for i=1:n2-1
        dx2_dt(i) = (x2(i+1) - x2(i)) / dt2(i);
        dy2_dt(i) = (y2(i+1) - y2(i)) / dt2(i);
    end
    
    dx1_dt(end,1) = dx1_dt(end-1,1);
    dy1_dt(end,1) = dy1_dt(end-1,1);
    dx2_dt(end,1) = dx2_dt(end-1,1);
    dy2_dt(end,1) = dy2_dt(end-1,1);
    
    dx1_noise = sigma_noise * randn(n,1);
    dy1_noise = sigma_noise * randn(n,1);
    
    dx2_noise = sigma_noise * randn(n,1);
    dy2_noise = sigma_noise * randn(n,1);
    
    dx1_dt = dx1_dt + dx1_noise;
    dy1_dt = dy1_dt + dy1_noise;
    
    x1(2:n) = x1(1:n-1) + dx1_dt(1:n-1) .* dt1(1:n-1);
    y1(2:n) = y1(1:n-1) + dy1_dt(1:n-1) .* dt1(1:n-1);
    
    dx2_dt = dx2_dt + dx2_noise;
    dy2_dt = dy2_dt + dy2_noise;
    
    x2(2:n) = x2(1:n-1) + dx2_dt(1:n-1) .* dt2(1:n-1);
    y2(2:n) = y2(1:n-1) + dy2_dt(1:n-1) .* dt2(1:n-1);
    
    tim1 = cumsum(dt1);
    tim2 = cumsum(dt2);
    
    tim1 = [0; tim1];
    tim2 = [0; tim2];
    
    % tau1 and tau2(partly)
    for i=1:n1-1
        x1_curr = x1(i); y1_curr = y1(i);
        vx1_curr = dx1_dt(i); vy1_curr = dy1_dt(i);
        
        ind2 = find(tim2 >= tim1(i)); %handle case when this is
                                      %empty
        if isempty(ind2)
            break;
        end
        
        j = ind2(1);
        x2_curr = x2(j); y2_curr = y2(j);
        vx2_curr = dx2_dt(j); vy2_curr = dy2_dt(j);
        
        tau_curr = findCollisionTime(x1_curr, y1_curr, vx1_curr, ...
                                     vy1_curr, x2_curr, y2_curr, ...
                                     vx2_curr, vy2_curr);
        
        tau1(i) = tau_curr; tau2(j) = tau_curr;
        
    end
    tau1(end) = 100;
    tau2(end) = 100;
    
    
    
end