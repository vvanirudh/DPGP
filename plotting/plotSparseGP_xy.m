function plotSparseGP_xy( x, y, dx_dt, dy_dt, data_x, data_y, sparseGP_x, sparseGP_y, act_reg)

global act_reg_g
dx_dt_zero = zeros(size(x));
dy_dt_zero = zeros(size(y));

x_min = min(x)-0.5;
x_max = max(x)+0.5;
y_min = min(y)-0.5;
y_max = max(y)+0.5;

% for plotting
x_range = x_max - x_min; 
y_range = y_max - y_min;
if x_range > y_range
   y_min = (y_min+y_max)/2.0 - x_range/2;
   y_max = (y_min+y_max)/2.0 + x_range/2;
else
   x_min = (x_min+x_max)/2.0 - y_range/2;
   x_max = (x_min+x_max)/2.0 + y_range/2;
end


% x velocity field
figure
%subplot(3,2,1)
%quiver(x,y,dx_dt,dy_dt_zero,'r');
%hold on
%plot (data_x,data_y, 'ro',sparseGP_x.BV_x, sparseGP_x.BV_y,'k*');
%axis([x_min x_max y_min y_max]);
%title('Sparse GP: x velocity field');
%axis equal
% y velcity field
%subplot(3,2,3)
%quiver(x,y,dx_dt_zero,dy_dt,'b');
%hold on
%plot (data_x,data_y, 'ro',sparseGP_y.BV_x, sparseGP_y.BV_y,'k*');
%axis([x_min x_max y_min y_max]);
%title('Sparse GP: y velocity field');
%axis equal
% x-y velocity field combined

%subplot(3,2,5)
quiver(x,y,dx_dt,dy_dt);
hold on
plot (data_x,data_y, 'ro');
axis([x_min x_max y_min y_max]);
title('Sparse GP: x-y velocity field');
axis equal

% variance

% calculate mean and variance
%x_query = x_min: 0.5: x_max;
%y_query = y_min: 0.5: y_max;
%[X_query,Y_query] = meshgrid(x_query, y_query);

%[n m] = size(X_query);
%var_x = ones(n,m);
%var_y = ones(n,m);

%for i = 1:n
%    for j = 1:m
%        [mu_x, var_x(i,j)] = sparseGP_predict(sparseGP_x, X_query(i,j), Y_query(i,j));
%        [mu_y, var_y(i,j)] = sparseGP_predict(sparseGP_y, X_query(i,j), Y_query(i,j));
%    end
%end
% sparseGP_x variance
%subplot(3,2,2)
%contour(X_query, Y_query, var_x,'ShowText','on');
%hold on
%plot (data_x,data_y, 'ro');
%axis([x_min x_max y_min y_max]);
%title('Sparse GP x: variance');
%axis equal
% sparseGP_y variance
%subplot(3,2,4)
%contour(X_query, Y_query, var_y,'ShowText','on');
%hold on
%plot (data_x,data_y, 'ro');
%axis([x_min x_max y_min y_max]);
%title('Sparse GP y: variance');
%axis equal
%set(gcf, 'Position', [1100, 100, 800, 800]);

% sparseGP_loc variance
%if (nargin == 9)

%    subplot(3,2,6)
%    [X_query,Y_query] = meshgrid(act_reg.x_ticks, act_reg.y_ticks);
%    prob_active = act_reg.obs_count;
%    [m, n] = size(prob_active);
%    for i = 1:m
%        for j = 1:n
%            if prob_active(i,j) > 0
%                count = prob_active(i,j);
%            else
%                count = act_reg.num_traj + 1;
%            end
%            act_reg_g.prob_active(count, act_reg.num_traj);
%        end
%    end
            
%    contour(X_query,Y_query,prob_active','ShowText','on');
%    hold on
%    %plot (data_x,data_y, 'ro');
%    axis([x_min x_max y_min y_max]);
%    title(sprintf('active region prediction with %d trajs', act_reg.num_traj));
%    %axis equal
%end



end