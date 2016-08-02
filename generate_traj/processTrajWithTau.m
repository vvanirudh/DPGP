function [t1, t2] = processTrajWithTau(x1, y1, dt1, x2, y2, dt2, ...
                                       sigma_noise)
   
    t1 = struct;
    t2 = struct;
    
    % Lengths of each trajectory
    n1 = length(x1);
    n2 = length(x2);
    
    % Initialize velocity and tau arrays
    dx1_dt = zeros(n1, 1);    
    dy1_dt = zeros(n1, 1);
    tau1 = zeros(n1, 1);
    
    dx2_dt = zeros(n2, 1);
    dy2_dt = zeros(n2, 1);
    tau2 = zeros(n2, 1);
    
    % Calculate velocities based on positions
    for i=1:n1-1
        dx1_dt(i) = (x1(i+1) - x1(i)) / dt1(i);
        dy1_dt(i) = (y1(i+1) - y1(i)) / dt1(i);
    end
    
    for i=1:n2-1
        dx2_dt(i) = (x2(i+1) - x2(i)) / dt2(i);
        dy2_dt(i) = (y2(i+1) - y2(i)) / dt2(i);
    end
    
    % Final velocity is equal to penultimate velocity
    dx1_dt(end,1) = dx1_dt(end-1,1);
    dy1_dt(end,1) = dy1_dt(end-1,1);
    dx2_dt(end,1) = dx2_dt(end-1,1);
    dy2_dt(end,1) = dy2_dt(end-1,1);
       
    
    tim1 = cumsum(dt1);
    tim2 = cumsum(dt2);
    
    tim1 = [0; tim1];
    tim2 = [0; tim2];
    
    % tau1 
    for i=1:n1-1
        x1_curr = x1(i); y1_curr = y1(i);
        vx1_curr = dx1_dt(i); vy1_curr = dy1_dt(i);
        
        ind2 = find(tim2(1:n2-1) <= tim1(i)); %handle case when this is
                                      %empty
        if isempty(ind2)
            break;
        end
        
        j = ind2(end);
        x2_curr = x2(j); y2_curr = y2(j);
        vx2_curr = dx2_dt(j); vy2_curr = dy2_dt(j);
        
        x2_curr = x2_curr + vx2_curr*(tim1(i) - tim2(j));
        y2_curr = y2_curr + vy2_curr*(tim1(i) - tim2(j));
        
        tau_curr = findCollisionTime(x1_curr, y1_curr, vx1_curr, ...
                                     vy1_curr, x2_curr, y2_curr, ...
                                     vx2_curr, vy2_curr);
        
        tau1(i) = tau_curr;
    end
    
    % tau2
    for i=1:n2-1
        x2_curr = x2(i); y2_curr = y2(i);
        vx2_curr = dx2_dt(i); vy2_curr = dy2_dt(i);
        
        ind1 = find(tim1(1:n1-1) <= tim2(i)); %handle case when this is
                                      %empty
        if isempty(ind1)
            break;
        end
        
        j = ind1(end);
        x1_curr = x1(j); y1_curr = y1(j);
        vx1_curr = dx1_dt(j); vy1_curr = dy1_dt(j);
        
        x1_curr = x1_curr + vx1_curr*(tim2(i) - tim1(j));
        y1_curr = y1_curr + vy1_curr*(tim2(i) - tim1(j));
        
        tau_curr = findCollisionTime(x2_curr, y2_curr, vx2_curr, ...
                                     vy2_curr, x1_curr, y1_curr, ...
                                     vx1_curr, vy1_curr);
        
        tau2(i) = tau_curr;
    end
    
    tau1(end) = 10;
    tau2(end) = 10;
    
    
    % Noise
    dx1_noise = sigma_noise * randn(n1,1);
    dy1_noise = sigma_noise * randn(n1,1);
    
    dx2_noise = sigma_noise * randn(n2,1);
    dy2_noise = sigma_noise * randn(n2,1);
    
    % Add noise to the velocities
    dx1_dt = dx1_dt + dx1_noise;
    dy1_dt = dy1_dt + dy1_noise;
    
    % Update the positions according to the noise
    x1(2:n1) = x1(1:n1-1) + dx1_dt(1:n1-1) .* dt1(1:n1-1);
    y1(2:n1) = y1(1:n1-1) + dy1_dt(1:n1-1) .* dt1(1:n1-1);
    
    % Add noise to the velocities
    dx2_dt = dx2_dt + dx2_noise;
    dy2_dt = dy2_dt + dy2_noise;
    
    % Update the positions according to the noise
    x2(2:n2) = x2(1:n2-1) + dx2_dt(1:n2-1) .* dt2(1:n2-1);
    y2(2:n2) = y2(1:n2-1) + dy2_dt(1:n2-1) .* dt2(1:n2-1);
    
    % Fill up the structures
    t1.x = x1; t1.y = y1; t1.dx_dt = dx1_dt; t1.dy_dt = dy1_dt;
    t1.tau = tau1;
    
    t2.x = x2; t2.y = y2; t2.dx_dt = dx2_dt; t2.dy_dt = dy2_dt;
    t2.tau = tau2;
    
end