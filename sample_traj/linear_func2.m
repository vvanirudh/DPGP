function [x,y] = linear_func2(n, pLimits)
    
    offset = 1;
    x_min = pLimits(1)+offset;
    x_max = pLimits(2)-offset;
    y_min = pLimits(3)+offset;
    y_max = pLimits(4)-offset;

    a = -0.5;
    b = 0.2;
    
    dx = (x_max - x_min) / n;
    dx_old = (1-(-1))/n;
    scale = dx/dx_old;
    x = (-1:dx_old:1)';
    y = a * x - b;


    x = x*scale;
    y = y*scale;
    
    x = flipud(x);
    y = flipud(y);
    
end