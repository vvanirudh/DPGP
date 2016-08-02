function plotSparseGP_xytau(x, y, tau, dx_dt, dy_dt, data_x, data_y, ...
                            data_tau, sparseGP_x, sparseGP_y, ...
                            act_reg)


    global act_reg_g
    
    dx_dt_zero = zeros(size(x));
    dy_dt_zero = zeros(size(y));
    dtau_dt_zero = zeros(size(tau));

    x_min = min(x)-0.5;
    x_max = max(x)+0.5;
    y_min = min(y)-0.5;
    y_max = max(y)+0.5;
    tau_min = min(tau) - 0.5;
    tau_max = max(tau) + 0.5;
    
    % for plotting
    x_range = x_max - x_min; 
    y_range = y_max - y_min;
    tau_range = tau_max - tau_min;
    if x_range > y_range
        y_min = (y_min+y_max)/2.0 - x_range/2;
        y_max = (y_min+y_max)/2.0 + x_range/2;
    else
        x_min = (x_min+x_max)/2.0 - y_range/2;
        x_max = (x_min+x_max)/2.0 + y_range/2;
    end
    
    figure
    quiver3(x, y, tau, dx_dt, dy_dt, dtau_dt_zero);
    hold on
    plot3(data_x,data_y, data_tau, 'ro');
    axis([x_min x_max y_min y_max tau_min tau_max]);
    title('Sparse GP: x-y-tau velocity field');
    axis equal


end
