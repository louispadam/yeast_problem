%% Boilerplate

% Add project to path
addpath(genpath(pwd));

% Build parameter structures
parameters = instantiate_parameters();

disp('Ran Boilerplate')

%% Set up problem

parameters.s1     = 0.00;
parameters.s2     = 0.10;
parameters.r1     = 0.20;
parameters.r2     = 0.30;
parameters.eps    = 0;
parameters.alph   = 0.9;
parameters.ct     = ct_sharp(0.1);

disp('Set Parameters')

%% Generate Initial Data

% Setup domain
x = linspace(0,1-2^(-16),2^16);

k = 4;
ic_s = @(z) (1-cos(2*pi*k*z))+0.5;
ic_t = ic_s(x);
s_n = sum([ic_t(2:end) + ic_t(1:end-1),ic_t(1) + ic_t(end)] .* ...
          [x(2:end) - x(1:end-1),mod(x(1) - x(end),1)])/2;
ic = @(z) ic_s(z)/s_n;

generate_sample = construct_sampler(ic);
ic_c = ic(x);

disp("Generated Initial Data")

%% Generate Parameter Space

param_space = generate_parameter_space((1:7)/10,(1:7)/10,(1:7)/10);

disp("Generated Parameter Space")

%% Simulate MF

fprintf("\nBegan Simulating MF\n");

% Experiment ID
id = "mf_data_26_10_15";

% Set up record of metadata
meta = fopen(id+'_metadata.txt', 'w');
fprintf(meta,"EXPERIMENT METADATA for ID: "+id+"\n"+...
             "===================\n\n");
fprintf(meta,"Date: %s\n\n",datetime('now'));
fprintf(meta,"Summary: Simulation results for MF across a wide "+...
             "parameter space. The goal is to collect MF data "+...
             "separately from NODE. In principle, MF need only be "+...
             "simulated once, and the results repeatedly referred "+...
             "to in the future. Simulation of NODE is a more "+...
             "fragile because of the sampling aspect, but also much "+...
             "faster and involves debugging a separate file.\n\n");
fprintf(meta,"System Parameters: Simulations were run until "+...
             "convergence (no maximum walltime), tracing 7 cluster "+...
             "modes using quantify_clustering and analyze_clustering. "+...
             "Initial conditions used an equipartitioned mesh of "+...
             "2^16 points on [0,1) with the density (after "+...
             "normalization) function f(x)=1.5-cos(8pi x). The space "+...
             "of model parameters was constructed using "+...
             "generate_parameter_space with uniform inputs (1:7)/10.");

fclose(meta);

% Run experiment
fprintf('Current time: %s\n', char(datetime('now','Format','HH:mm:ss')));
exp_data = collect_mf_data(parameters,ic_c,param_space,...
                           'Modes',7,...
                           'Max_Simulation_Time',Inf,...
                           'Save_Data',id);

fprintf("Finished Big Experiment\n");

%% Run Big Experiment

fprintf("\nBegan Big Experiment\n");

fprintf('Current time: %s\n', char(datetime('now','Format','HH:mm:ss')));

exp_data = big_experiment(parameters,ic_c,generate_sample,param_space,...
                          'Modes',7,...
                          'Max_Simulation_Time',2000,...
                          'Trials',10,...
                          'N',1000:1000:5000,...
                          'Save_Data',"big_experiment_26_10_06");

fprintf("Finished Big Experiment\n");

