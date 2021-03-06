%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementation of 
% DPGP (A Baysesian Nonparametric Approach to 
%   Modeling Motion Patterns) 
%   by Joseph, Roy, et al.
% Steven Chen, Dec 21, 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;
close all;
clc;

% adding directory path
addpath ..
addpath ../plotting
addpath ../generate_traj
addpath ../sample_traj
addpath ../data
addpath ../Gibbs
% global variables
global trajs sparseGPs
lx = 2;
ly = 2;


%% generate n trajectories
% to do: implement generateTrajs
n_traj = 50;
n_points = 30;
x_min = -5; x_max = 5; y_min = -5; y_max = 5;
pLimit = [x_min, x_max, y_min, y_max];
speed = 1.0;
sigma_noise = 0.05;
% rng(1);
trajs = generateTrajsVaryingSpeed(n_traj, n_points, pLimit, speed, sigma_noise);

n_traj = trajs.n_traj;
plotTrajs(trajs, 'initialization');

%% DPGP
sigma_noise = 1.0;
sigma_input = 1;
hyperparam = [lx, ly, sigma_input, sigma_noise];
n_sweep = 200;
% randomly assign trajectories to clusters
% to do: implement initialization routine
% to do: keep track of all samples (Gibbs sampling)
trajs.n_clus = round(log(n_traj));
cluster = zeros(n_traj, n_sweep);
cluster(:,1) = unidrnd(trajs.n_clus, n_traj,1);
%for debugging, start with perfect initialization
%for i = 1:n_traj
%   cluster(i,1) = mod(i,4)+1; 
%end

trajs.cluster = cluster;
count = groupTraj(1);
alpha = 0.5;
initialize_SparseGPs_array(hyperparam);
%count = groupTraj(1);
%build_SparseGPs_array(hyperparam);
%plotTrajs(trajs);
%plotSparseGP_array(sparseGPs,5);

%% main loop
tic
for sweep_num = 1:n_sweep
  trajs.sweep_count = sweep_num;
  % construct sparseGP for each motion pattern
  % to do: call sparseGP constructor
  if (sweep_num > 1)
    build_SparseGPs_array(hyperparam);
  end
  % for initialization
  %if (sweep_num == 1)
  %    trajs.cluster(trajs.n_clus+1:end,1) = abs(trajs.cluster(trajs.n_clus+1:end,1));
  %end
  

  % update hyperparam \alpha
  % to do: implement updateAlpha
  % alpha = update_alpha(trajs.n_clus, trajs.n_traj);
    fprintf(sprintf('%d th sweep, %d clusters, alpha: %.2f \n',sweep_num, trajs.n_clus, alpha));
  for pp =1:length(count)
    fprintf(sprintf('%d ',count(pp)));
  end
  fprintf(sprintf('\n'));
  
  
  % for each trajectories
  % compute likelihood of drawing from each existing cluster
  L = zeros(trajs.n_traj, trajs.n_clus+1);
  L_GP = zeros(trajs.n_traj, trajs.n_clus+1);
  rand_ordering = 1:trajs.n_traj; randperm(trajs.n_traj);

  for kk = 1:trajs.n_traj % kth traj
      k = rand_ordering(kk);
      % existing clusters
      for j = 1:trajs.n_clus % jth cluster
          trajs.data(k).DP_alpha = alpha;
          if j == trajs.cluster(k,sweep_num);
              n_k = count(j)-1;
          else
              n_k = count(j);
          end

          if (n_k == 0)
              n_k = alpha;
          end
          
          %else
              %L_GP(k, j) = DP_traj_likelihood(sparseGPs(j).sparseGP_x, sparseGPs(j).sparseGP_y, trajs.data(k));
          L_GP(k, j) = DP_traj_likelihood_indep(sparseGPs(j).sparseGP_x, sparseGPs(j).sparseGP_y, trajs.data(k));
          L(k, j) = L_GP(k, j) + log (n_k / (trajs.n_traj-1+alpha));
          %end
          %fprintf(sprintf('sweep: %d, k: %d, j: %d, n_k: %d \n', sweep_num, k ,j, n_k));
      end

      
      normalization_const = max(L(k,1:end-1));
      l = L(k,:)- normalization_const; % normalization
      l_GP = L_GP(k,:)- normalization_const; % normalization

      p = exp(l);
      p(trajs.n_clus+1) = mean(exp(l_GP(1:end-1))) * alpha / (trajs.n_traj-1+alpha);

      z_k = randsample(trajs.n_clus+1, 1, true, p);
      % update cluster
      if (z_k > trajs.n_clus)
          z_k = length(count) + 1; % new cluster
          count = [count; 0];
      end
      count(trajs.cluster(k,sweep_num)) = count(trajs.cluster(k,sweep_num)) - 1;
      trajs.cluster(k,sweep_num) = z_k;
      count(z_k) = count(z_k) + 1;
  end
  % initialize for the next round
  count = groupTraj(sweep_num);
  if sweep_num < n_sweep
      trajs.cluster(:,sweep_num+1) = trajs.cluster(:,sweep_num);
  end
  
    
    
end
timeElapsed = toc;
fprintf(sprintf('time elasped: %fs',timeElapsed));

%plotTrajs(trajs);
%% post processing
% burning in
% splicing


% calculate mean
%% plotting
burn_in = floor(n_sweep / 2);
splicing = 2;
[avgSample, mode, config, config_count] = gibbs_sampling_postProcessing(trajs.cluster', burn_in, splicing);
trajs.cluster(:,end) = mode';
plotTrajs(trajs, 'mode');

trajs.cluster(:,end) = avgSample';
plotTrajs(trajs, 'average sample');

% reassgined trajectories based on their current cluster config
groupTraj(sweep_num);
build_SparseGPs_array(hyperparam);
for i = 1:trajs.n_traj
    ind_order = findBestPattern(trajs.data(i));
   trajs.cluster(i,end) = ind_order(1);
end
plotTrajs(trajs, 'reassigned cluster');

build_SparseGPs_array(hyperparam);


count = groupTraj(sweep_num);
build_SparseGPs_array(hyperparam);
plotSparseGP_array(sparseGPs,5);
mode

% history of DP samples
figure
plot(config_count, 'r-o');
xlabel('DP sampled configuration label');
ylabel('counts')
