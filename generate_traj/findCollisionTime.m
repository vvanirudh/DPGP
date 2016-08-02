function tau = findCollisionTime(x1, y1, dx1_dt, dy1_dt, x2, y2, ...
                                 dx2_dt, dy2_dt)

    t = 0;
    tstep = 0.1;
    dist = Inf;
    
    while dist > 0.5 && t < 4.99
       
        xc1 = x1 + dx1_dt * t;
        yc1 = y1 + dy1_dt * t;
        
        xc2 = x2 + dx2_dt * t;
        yc2 = y2 + dy2_dt * t;
        
        dist = min(dist, norm([xc1, yc1] - [xc2, yc2]));
        
        if dist < 0.5
            break;
        end
        
        t = t + tstep;
    end
    
    tau = t;
end
