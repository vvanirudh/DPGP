function [x, y] = cubic(n, pLimits)

    offset = 1;
    x_min = pLimits(1)+offset;
    x_max = pLimits(2)-offset;
    y_min = pLimits(3)+offset;
    y_max = pLimits(4)-offset;

    dx = (x_max - x_min) / n;
    dx_old = (1-(-1))/n;
    scale = dx/dx_old;
    x = (-1:dx_old:1)';
    y = x.^3;


    x = x*scale;
    y = y*scale;

end