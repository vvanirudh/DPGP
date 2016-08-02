function [ind_order, ifNew] = findBestPattern_tau(traj)

    global sparseGPs
    n = numel(sparseGPs);


    L_GP = zeros(1, n);
    L = zeros(1, n);
    num_counts = zeros(1,n);
    alpha = 0;
    ifNew = 0;

    % find the best fitting cluster
    for i = 1:n
        % sparse GP
        L_GP(i) = DP_traj_likelihood_indep_tau(sparseGPs(i).sparseGP_x, sparseGPs(i).sparseGP_y, traj);
        % active region
        % L_GP(i)= activeRegion_likelihoodS(sparseGPs(i).act_reg, traj);
        
        % DP component
        num_counts(i) = sparseGPs(i).count;
        alpha = traj.DP_alpha;  
    end

    % L_GP = L_GP;
    L = L_GP + log(num_counts / (sum(num_counts)+alpha))/3;
    
    [~, ind_order] = sort(L, 'descend');
    
    for k = 1:n
        if (L_GP(ind_order(k))/length(traj.x)) < -1.0%log(0.3)
            ifNew = k;
            break;
        end
    end



end
