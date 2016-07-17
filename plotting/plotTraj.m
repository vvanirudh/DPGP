function plotTraj( x, y, dx_dt, dy_dt )

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
plot (x,y, 'ro');
axis([x_min x_max y_min y_max]);
title('raw data: x velocity field');

% y velcity field
subplot(3,1,2)
quiver(x,y,dx_dt_zero,dy_dt);
hold on
plot (x,y, 'ro');
axis([x_min x_max y_min y_max]);
title('raw data: y velocity field');

% x-y velocity field combind
subplot(3,1,3)
quiver(x,y,dx_dt,dy_dt);
hold on
plot (x,y, 'ro');
axis([x_min x_max y_min y_max]);
title('raw data: x-y velocity field');

set(gcf, 'Position', [100, 100, 400, 800]);

end

