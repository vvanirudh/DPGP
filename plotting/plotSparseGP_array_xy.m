function plotSparseGP_array_xy(sparseGPs, max_shown_pattern, ifLoc_in)

if nargin == 3
    ifLoc = ifLoc_in;
else 
    ifLoc = 0;
end

n = min(max_shown_pattern, length(sparseGPs));
for i = 1:n
    x_min = min(sparseGPs(i).data.x);
    x_max = max(sparseGPs(i).data.x);
    y_min = min(sparseGPs(i).data.y);
    y_max = max(sparseGPs(i).data.y);
    
    x_query = x_min: 0.5: x_max;
    y_query = y_min: 0.5: y_max;
    [X_query,Y_query] = meshgrid(x_query, y_query);
    X_query = reshape(X_query,[length(x_query)*length(y_query),1]);
    Y_query = reshape(Y_query,[length(x_query)*length(y_query),1]);
    [x_vel_sparse, x_var_sparse] = sparseGP_predict(sparseGPs(i).sparseGP_x, X_query, Y_query);
    [y_vel_sparse, y_var_sparse] = sparseGP_predict(sparseGPs(i).sparseGP_y, X_query, Y_query);
    if (ifLoc == 0)
        plotSparseGP_xy( X_query, Y_query, x_vel_sparse, y_vel_sparse,  sparseGPs(i).data.x, sparseGPs(i).data.y, ...
            sparseGPs(i).sparseGP_x, sparseGPs(i).sparseGP_y)
    else
        plotSparseGP_xy( X_query, Y_query, x_vel_sparse, y_vel_sparse,  sparseGPs(i).data.x, sparseGPs(i).data.y, ...
            sparseGPs(i).sparseGP_x, sparseGPs(i).sparseGP_y, sparseGPs(i).act_reg)
    end
    
end

end


