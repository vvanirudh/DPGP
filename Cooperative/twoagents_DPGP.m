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

% global variables
global trajs sparseGPs

% hyperparameters (some of them)
lx = 2;
ly = 2;

%% generate n trajectories for both agents
n_traj = 50;
n_points = 30;
x_min = -5; x_max = 5; y_min = -5; y_max = 5;
pLimit = [x_min, x_max, y_min, y_max];
speed = 1.0;
sigma_noise = 0.05;

% Implement generateTwoAgentTrajs
trajs = generateTwoAgentTrajs(n_traj, n_points, pLimit, speed, ...
                              sigma_noise);

n_traj = trajs.n_traj;

%% Two Agent DPGP
sigma_noise = 1.0;
sigma_input = 1;
hyperparam = [lx, ly, sigma_input, sigma_noise];
n_sweep = 200;

% Randomly assign trajectories to clusters
trajs.n_clus = round(log(n_traj))*2; % * 2 because of two agent
                                     % trajectories
cluster = zeros(n_traj, 2, n_sweep);
cluster(:, 1) = unidrnd(trajs.n_clus, n_traj, 2, 1);

trajs.cluster = cluster;
