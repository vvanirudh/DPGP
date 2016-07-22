%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trying out to come up with a cooperative
% modeling for two agents using the 
% DPGP Model
% by Anirudh Vemula, Jul 18, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

% adding directory path
addpath ..
addpath ../plotting
addpath ../generate_traj
addpath ../sample_traj
addpath ../Gibbs
addpath ../DP
addpath genTraj
addpath ../external
addpath twoagents_plotting

% global variables
global trajs sparseGPs

% hyperparameters (some of them)
lx = 2;
ly = 2;

%% generate n trajectories for both agents
n_traj = 7;
n_points = 30;
x_min = -5; x_max = 5; y_min = -5; y_max = 5;
pLimit = [x_min, x_max, y_min, y_max];
speed = 1.0;
sigma_noise = 0.05;

% Implement generateTwoAgentTrajs
trajs = generateTwoAgentTrajs(n_traj, n_points, pLimit, speed, ...
                              sigma_noise);

n_traj = trajs.n_traj;
%plotTrajs(trajs, 'initialization');

%% Two Agent DPGP
sigma_noise = 1.0;
sigma_input = 1;
hyperparam = [lx, ly, sigma_input, sigma_noise];
n_sweep = 200;

% Randomly assign trajectories to clusters
trajs.n_clus = round(log(n_traj))*2; % * 2 because of two agent
                                     % trajectories
cluster = zeros(n_traj, 2, n_sweep);
cluster(:, :, 1) = unidrnd(trajs.n_clus, n_traj, 2, 1);

trajs.cluster = cluster;
count = twoagents_groupTraj(1);
alpha = 0.5;

% Initialize and fit the GPs to the current data
twoagents_initialize_SparseGPs_array(hyperparam);

%% main loop
tic
for sweep_num=1:n_sweep
    trajs.sweep_count = sweep_num;
    % Construct sparseGP for each motion pattern
    if (sweep_num > 1)
        twoagents_build_SparseGPs_array(hyperparam);
    end
    
    fprintf(sprintf('%d th sweep, %d clusters, alpha: %.2f\n', ...
                    sweep_num, trajs.n_clus, alpha));
    for pp=1:length(count)
        fprintf(sprintf('%d ', count(pp)));
    end
    fprintf(sprintf('\n'));
    
    % For each pair of trajectories
    % Compute *JOINT* likelihood of drawing from each existing cluster
    
    L = zeros(trajs.n_traj, trajs.n_clus+1, trajs.n_clus+1);
    L_GP = zeros(trajs.n_traj, trajs.n_clus+1, trajs.n_clus+1);
    rand_ordering = 1:trajs.n_traj; %randperm(trajs.n_traj)
    
    for kk=1:trajs.n_traj
        k = rand_ordering(kk);
        
        % p(t1, t2) = p(t1)*p(t2|t1)
        % Initially, likelihood p(t|b)
        for j1=1:trajs.n_clus
            trajs.data1(k).DP_alpha = alpha;
            
            if j1 == trajs.cluster(k, 1, sweep_num)
                n_k_1 = count(j1) - 1;
            else
                n_k_1 = count(j1);
            end
            
            if(n_k_1 == 0)
                n_k_1 = alpha;
            end
            
            % Calculate the likelihood of p(t1|b_j1)
            likelihood_1_GP = ...
                DP_traj_likelihood_indep(sparseGPs(j1).sparseGP_x, ...
                                         sparseGPs(j1).sparseGP_y, ...
                                         trajs.data1(k));
            likelihood_1 = likelihood_1_GP + log (n_k_1 / ...
                                                  (trajs.n_traj*2 ...
                                                   - 1 + alpha));
            
            for j2=1:trajs.n_clus                
                trajs.data2(k).DP_alpha = alpha;
                
                if j2 == trajs.cluster(k, 2, sweep_num)
                    n_k_2 = count(j2) - 1;
                else
                    n_k_2 = count(j2);
                end                                
                
                if(n_k_2 == 0)
                    n_k_2 = alpha;
                end
                
                % Now to calculate the likelihood of p(t2|b_j1, t1,
                % b_j2)
                likelihood_2_GP = ...
                    DP_traj_likelihood_dep(sparseGPs(j2).sparseGP_x, ...
                                           sparseGPs(j2).sparseGP_y, ...
                                           trajs.data2(k), ...
                                           trajs.data1(k));                
                
                likelihood_2 = likelihood_2_GP + log(n_k_2 / ...
                                                     (trajs.n_traj*2 ...
                                                      - 1 + ...
                                                      alpha));
                
                L_GP(k, j1, j2) = likelihood_1_GP + ...
                    likelihood_2_GP;
                L(k, j1, j2) = likelihood_1 + likelihood_2;                    
            end
        end
        
        % DOUBT
        normalization_const = max(max(L(k, 1:end-1, 1:end-1)));
        l = squeeze(L(k,:,:)) - normalization_const;
        l_GP = squeeze(L_GP(k,:,:)) - normalization_const;

        p = exp(l);
        p(trajs.n_clus+1, 1:end-1) = mean(exp(l_GP(1:end-1, 1:end-1))) ...
            * alpha / (trajs.n_traj*2 - 1 + alpha);
        
        p(1:end-1, trajs.n_clus + 1) = mean(exp(l_GP(1:end-1, 1: ...
                                                     end-1)), 2) * ...
            alpha / (trajs.n_traj*2 - 1 + alpha);
        
        p(trajs.n_clus+1, trajs.n_clus+1) = mean2(exp(l_GP(1:end-1, ...
                                                          1:end-1))) ...
            * (alpha / (trajs.n_traj*2 - 1 + alpha))^2;
        
        %keyboard()
        
        % Draw z_k_1 and z_k_2 from the probability distribution
        [z_k_1 z_k_2] = pinky(1:trajs.n_clus+1, 1:trajs.n_clus+1, ...
                              p');
        
        % If both are assigned to the new cluster
        if(z_k_1 > trajs.n_clus && z_k_2 > trajs.n_clus)
            z_k_1 = length(count) + 1;
            z_k_2 = length(count) + 1;
            count = [count ; 0];           
        elseif (z_k_1 > trajs.n_clus) 
            % Only the 1st one
            z_k_1 = length(count) + 1;
            count = [count ; 0];
        elseif (z_k_2 > trajs.n_clus)
            % Only the second one
            z_k_2 = length(count) + 1;
            count = [count ; 0];
        end
        
        count(trajs.cluster(k, 1, sweep_num)) = count(trajs.cluster(k, ...
                                                          1, sweep_num)) ...
            - 1;
        count(trajs.cluster(k, 2, sweep_num)) = count(trajs.cluster(k, ...
                                                          2, sweep_num)) ...
            - 1;
        
        trajs.cluster(k, 1, sweep_num) = z_k_1;
        trajs.cluster(k, 2, sweep_num) = z_k_2;
        count(z_k_1) = count(z_k_1) + 1;
        count(z_k_2) = count(z_k_2) + 1;
        
    end
    
    %Initialize for the next round
    count = twoagents_groupTraj(sweep_num);
    if sweep_num < n_sweep
        trajs.cluster(:,:,sweep_num+1) = trajs.cluster(:, :, sweep_num);
    end
end

timeElapsed = toc;
fprintf(sprintf('time elasped: %fs',timeElapsed));


% calculate mean
%% plotting

burn_in = floor(n_sweep / 2);
splicing = 2;
[avgSample1, mode1, config1, config_count1] = ...
    gibbs_sampling_postProcessing(squeeze(trajs.cluster(:,1,:))', ...
                                  burn_in, splicing);
[avgSample2, mode2, config2, config_count2] = ...
    gibbs_sampling_postProcessing(squeeze(trajs.cluster(:,2,:))', ...
                                  burn_in, splicing);

trajs.cluster(:, 1, end) = mode1';
trajs.cluster(:, 2, end) = mode2';

%keyboard()
twoagents_plotTrajs(trajs, 'mode');

trajs.cluster(:, 1, end) = avgSample1';
trajs.cluster(:, 2, end) = avgSample2';

twoagents_plotTrajs(trajs, 'avgSample');

%twoagents_groupTraj(sweep_num);
%twoagents_build_SparseGPs_array(hyperparam);
%for i=1:trajs.n_traj
%    ind_order = twoagents_findBestPattern(trajs.data1(i), ...
%                                          trajs.data2(i));
%    trajs.cluster(i, 1, end) = ind_order(1, 1);
%    trajs.cluster(i, 2, end) = ind_order(2, 1);
%end