%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to find the best cluster pair
% that fits the current pair of trajectories
% by Anirudh Vemula, Jul 21, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ind1, ind2] = twoagents_findBestPattern(traj1, ...
                                                  traj2, ifIndep)

    global sparseGPs
    n = numel(sparseGPs);
    
    L_GP = zeros(n, n);
    L = zeros(n, n);
    
    num_counts = zeros(n, n);
    alpha = 0;
    
    % Obtain the total number of trajectories
    numTrajs = 0;
    for i=1:n
        numTrajs = numTrajs + sparseGPs(i).count;
    end
    
    % Get alpha value
    alpha = traj1.DP_alpha;
    
    % find the best fitting pair of clusters
    for i=1:n
        for j=1:n
            % sparse GP
            if ~ifIndep
                L_GP(i, j) = DP_traj_likelihood_indep(sparseGPs(i).sparseGP_x, ...
                                                      sparseGPs(i).sparseGP_y, ...
                                                      traj1) + ...
                    DP_traj_likelihood_dep(sparseGPs(j).sparseGP_x, ...
                                           sparseGPs(j).sparseGP_y, ...
                                           traj2, traj1);
            else
                L_GP(i, j) = DP_traj_likelihood_indep(sparseGPs(i).sparseGP_x, ...
                                                      sparseGPs(i).sparseGP_y, ...
                                                      traj1) + ...
                    DP_traj_likelihood_indep(sparseGPs(j).sparseGP_x, ...
                                             sparseGPs(j).sparseGP_y, ...
                                             traj2);
            end
            
            num_counts(i, j) = log (sparseGPs(i).count / (numTrajs ...
                                                          + alpha))/3 ...
                + log (sparseGPs(j).count / (numTrajs + alpha))/3;                        
        end
    end
            
    L = L_GP + num_counts;
    
    [ind1, ind2] = ind2sub(size(L), find(L == max(max(L))));
    
end