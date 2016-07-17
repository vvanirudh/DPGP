function plotGP( x, y, dx_dt, dy_dt, data_x, data_y )

dx_dt_zero = zeros(size(x));
dy_dt_zero = zeros(size(y));

x_min = min(x)-1;
x_max = max(x)+1;
y_min = min(y)-1;
y_max = max(y)+1;

% x velocity field
figure
subplot(3,1,1)
quiver(x,y,dx_dt,dy_dt_zero,'r');
hold on
plot (data_x,data_y, 'ro')
axis([x_min x_max y_min y_max]);
title('Batch GP: x velocity field');

% y velcity field
subplot(3,1,2)
quiver(x,y,dx_dt_zero,dy_dt,'b');
hold on
plot (data_x,data_y, 'ro')
axis([x_min x_max y_min y_max]);
title('Batch GP: y velocity field');

% x-y velocity field combind
subplot(3,1,3)
quiver(x,y,dx_dt,dy_dt,1);
hold on
plot (data_x,data_y, 'ro')
axis([x_min x_max y_min y_max]);
title('Batch GP: x-y velocity field');

set(gcf, 'Position', [600, 100, 400, 800]);


% for debugging
for i=2:length(x)
    if (y(i) == y(1))
        n = i - 1;
        m = length(x)/n;
        break;
    end
end
%m, n
%reshape(dx_dt,[m,n])
% figure 
% surf(reshape(x,[m,n]),reshape(y,[m,n]),reshape(dx_dt,[m,n]));
% figure
% surf(reshape(x,[m,n]),reshape(y,[m,n]),reshape(dy_dt,[m,n]));

end