%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to find the best cluster pair
% that fits the current pair of trajectories
% by Anirudh Vemula, Jul 21, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ind_order, ifNew] = twoagents_findBestPattern(traj1, ...
                                                  traj2)

    global sparseGPs
    n = numel(sparseGPs);
    
    L_GP = zeros(n, n);
    L = zeros(n, n);
    
    num_counts = zeros(n, n);
    alpha = 0;
    ifNew = 0;
    
    % find the best fitting pair of clusters
    for i=1:n
        for j=1:n
            % sparse GP
            L_GP(i, j) = DP_traj_likelihood_indep(sparseGPs(i).sparseGP_x, ...
                                                  sparseGPs(i).sparseGP_y, ...
                                                  traj1) + ...
                DP_traj_likelihood_dep(sparseGPs(j).sparseGP_x, ...
                                       sparseGPs(j).sparseGP_y, ...
                                       traj2, traj1);
            
            num_counts(i,j)= sparseGPs()
    

end