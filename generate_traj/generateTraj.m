function [x, y, dt] = generateTraj(n, pLimit, speed, func)

    [x, y] = func(n, pLimit);
    
    % random offset
    x = x + (rand - 0.5);
    y = y + (rand - 0.5);

    dt = zeros(n+1,1);
    for i=1:n-1
        dt(i) = norm([x(i+1)-x(i) y(i+1)-y(i)]) / speed;
    end
    dt(n) = dt(n-1);

end
