function predPlot(pred_traj, sparseGP_x, sparseGP_y)

addpath ../

x = pred_traj(:,1);
y = pred_traj(:,2);
sigma_x = sqrt(pred_traj(:,3));
sigma_y = sqrt(pred_traj(:,4));

x_min = min(sparseGP_x.x_data)-1;
x_max = max(sparseGP_x.x_data)+1;
y_min = min(sparseGP_x.y_data)-1;
y_max = max(sparseGP_x.y_data)+1;

% updating function mean and variance
x_query = x_min: 0.5: x_max;
y_query = y_min: 0.5: y_max;
[X_query,Y_query] = meshgrid(x_query, y_query);

X_query = reshape(X_query,[length(x_query)*length(y_query),1]);
Y_query = reshape(Y_query,[length(x_query)*length(y_query),1]);
[dx_dt, var_x_sparse] = sparseGP_predict(sparseGP_x, X_query, Y_query);
[dy_dt, var_y_sparse] = sparseGP_predict(sparseGP_y, X_query, Y_query);
data_x = sparseGP_x.x_data;
data_y = sparseGP_x.y_data;

% x-y velocity field combined
figure
quiver(X_query,Y_query,dx_dt,dy_dt);
hold on
plot (data_x,data_y, 'ro');
%axis([x_min x_max y_min y_max]);
title('Predicted trajectory');
hold on
plot(x,y,'g-o');
legend('predictive flow field', 'training data', 'predicted traj'); 
% plot uncertainty ellipse
hold on 
pred_length = length(pred_traj(:,1));
ec = ellipseCord(pred_traj(1,1), pred_traj(1,2), pred_traj(1,3), pred_traj(1,4));
plot (ec(1,:), ec(2,:),'g');
last_point = pred_traj(1,1:2);
for i = 1:pred_length
    if norm(pred_traj(i,1:2) - last_point, 2) > 1.5
        ec = ellipseCord(pred_traj(i,1), pred_traj(i,2), pred_traj(i,3), pred_traj(i,4));
        last_point = pred_traj(i,1:2);
        plot (ec(1,:), ec(2,:),'g');
    end        
end
axis equal

end





function ec = ellipseCord(xc, yc, a,b)
x = linspace(-a,a,20);
y = b * sqrt(1 - x.^2 /a^2);
x = [x, flip(x), x(1)]+xc;
y = [y, -flip(y), y(1)]+yc;
ec = [x; y];
end