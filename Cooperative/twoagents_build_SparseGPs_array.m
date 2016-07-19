function twoagents_build_SparseGPs_array(hyperparam)

global sparseGPs

for i=1:length(sparseGPs)
   sparseGPs(i).sparseGP_x = build_sparseGP(sparseGPs(i).data.x, ...
                                            sparseGPs(i).data.y, ...
                                            sparseGPs(i).data.dx_dt, ...
                                            hyperparam);
   sparseGPs(i).sparseGP_y = build_sparseGP(sparseGPs(i).data.x, ...
                                            sparseGPs(i).data.y, ...
                                            sparseGPs(i).data.dy_dt, ...
                                            hyperparam);
end

%% debugging plot
% x_min = -5; x_max = 5; y_min = -5; y_max = 5;
% x_query = x_min: 0.5: x_max;
% y_query = y_min: 0.5: y_max;
% [X_query,Y_query] = meshgrid(x_query, y_query);
% 
% X_query = reshape(X_query,[length(x_query)*length(y_query),1]);
% Y_query = reshape(Y_query,[length(x_query)*length(y_query),1]);
% 
% [x_vel_sparse, x_var_sparse] = sparseGP_predict(sparseGPs(1).sparseGP_x, X_query, Y_query);
% [y_vel_sparse, y_var_sparse] = sparseGP_predict(sparseGPs(1).sparseGP_y, X_query, Y_query);
% plotSparseGP( X_query, Y_query, x_vel_sparse, y_vel_sparse, ...
%     sparseGPs(1).data.x, sparseGPs(1).data.y, ...
%     sparseGPs(1).sparseGP_x, sparseGPs(1).sparseGP_y);

end