%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementation of 
% Sparse Online Gaussian Processes
%       by Casto and Opper
%
% Steven Chen, Nov 3, 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;

addpath plotting
addpath sample_traj
addpath generate_traj
% include all global variables
% global_var;
% parameters for generating observations
n = 30;
x_min = -5; x_max = 5; y_min = -5; y_max = 5;
pLimit = [x_min, x_max, y_min, y_max];
speed = 0.5;

% set hyperparameter

lx = 5;
ly = 5;

sigma_input = 1;
sigma_noise = 0.05;
hyperparam = [lx, ly, sigma_input, sigma_noise];
sparseGP.hyperparam = hyperparam;

% generating observations1
rng(1)
[x, y, dt] = generateTraj(n, pLimit, speed, @quadratic);

%[x, y, dt] = generateTraj(n, pLimit, speed, @sinusoidal);
%[x, y, dt] = generateTraj(n, pLimit, speed, @linear_func);
%[x, y, dt] = generateTraj(n, pLimit, speed, @cubic);
[x, y, dx_dt, dy_dt] = processTraj(x, y, dt,  sigma_noise);
plot(x,y,'ro');
plotTraj( x, y, dx_dt, dy_dt );

% updating function mean and variance
x_query = x_min: 0.5: x_max;
y_query = y_min: 0.5: y_max;

[X_query,Y_query] = meshgrid(x_query, y_query);


X_query = reshape(X_query,[length(x_query)*length(y_query),1]);
Y_query = reshape(Y_query,[length(x_query)*length(y_query),1]);


% GP prediction
%hyperparam_x = optimize_GP_l(x, y, dx_dt, hyperparam);
%hyperparam_y = optimize_GP_l(x, y, dy_dt, hyperparam);
hyperparam_x = hyperparam;
hyperparam_y = hyperparam;
[x_vel, var_x] = GP_predict(x, y, dx_dt, X_query, Y_query, hyperparam_x);
[y_vel, var_y] = GP_predict(x, y, dy_dt, X_query, Y_query, hyperparam_y);
plotGP( X_query, Y_query, x_vel, y_vel, x, y);

%combine x and y directions
sparseGP_x = build_sparseGP(x, y, dx_dt, hyperparam);
sparseGP_y = build_sparseGP(x, y, dy_dt, hyperparam);
[x_vel_sparse, var_x_sparse] = sparseGP_predict(sparseGP_x, X_query, Y_query);
[y_vel_sparse, var_y_sparse] = sparseGP_predict(sparseGP_y, X_query, Y_query);
plotSparseGP( X_query, Y_query, x_vel_sparse, y_vel_sparse, x, y, sparseGP_x, sparseGP_y);

% debugging sparseGP_predict_distri.m
pos = [0,0];
pos_var = [1 0; 0 1];
[mu_dis, var_dis]= sparseGP_predict_distri(sparseGP_x, pos, pos_var);
[x_vel_sparse, var_x_sparse] = sparseGP_predict(sparseGP_x, pos(1), pos(2));

% predicting future trajectory given starting point
starting_point = [-2,0];
num_steps = 10;
dt = ones(num_steps,1);
pos_var_in = [0.5, 0; 0, 0.5];
pred_traj = sparseGP_trajPredict(sparseGP_x, sparseGP_y, starting_point, pos_var_in, dt, num_steps, speed);
predPlot(pred_traj, sparseGP_x, sparseGP_y);
