%% Boilerplate

% Add project to path
addpath(genpath(pwd));

% Build parameter structures
parameters = instantiate_parameters();

disp('Ran Boilerplate')

%% Set up problem

parameters.s1     = 0.00;%0.10;
parameters.s2     = 0.10;%0.30;
parameters.r1     = 0.55;%0.65;
parameters.r2     = 0.70;%0.80;
parameters.eps    = 0;
parameters.alph   = 0.9;
parameters.ct     = ct_sharp(0.01);

disp('Set Parameters')

%% Generate Initial Data

% Setup domain
x = linspace(0,1-2^(-16),2^16);

k = 4;
ic_func = @(z) (1-cos(2*pi*k*z))+0.5;
ic_eval = ic_func(x);
for_norm = sum([ic_eval(2:end) + ic_eval(1:end-1),...
                ic_eval(1) + ic_eval(end)] .* ...
               [x(2:end) - x(1:end-1),mod(x(1) - x(end),1)])/2;
ic_func_norm = @(z) ic_func(z)/for_norm;

generate_sample = construct_sampler(ic_func_norm);
ic_eval_norm = ic_func_norm(x);

<<<<<<< Updated upstream
ic_lag = euler_to_lagrange(x,ic_eval_norm);
ic_mass = ic_lag(2,:);

=======
>>>>>>> Stashed changes
disp('Prepared Initial Data')

%% Generate Parameter Space

base_set = [0,0.1,0.55,0.7];
param_space = generate_parameter_space_abs(base_set(1),...
                                           base_set(2) + (-1:1)/20,...
                                           base_set(3) + (-2:2)/20,...
                                           base_set(4) + (-2:2)/20);

disp("Generated Parameter Space")

%% MEFI EXPERIMENT

fprintf("\nCollecting MEFI Data\n");

% Experiment ID
id = "mefi_data_26_11_16";

% Set up record of metadata
meta = fopen(id+'_metadata.txt', 'w');
fprintf(meta,"EXPERIMENT METADATA for ID: "+id+"\n"+...
             "===================\n\n");
fprintf(meta,"Date: %s\n\n",datetime('now'));
fprintf(meta,"Summary: Simulation results for MEFI with parameter "+...
             "set [0,0.1,0.55,0.7], where we observe a transition "+...
             "from a 4-cluster state to a 2-cluster state.\n\n");
fprintf(meta,"System Parameters: Simulation was run until "+...
             "convergence (no maximum walltime), tracing 7 cluster "+...
             "modes. Initial conditions used an equipartitioned mesh "+...
             "of 2^16 points on [0,1) with the density (after "+...
             "normalization) function f(x)=1.5-cos(8pi x).");

fclose(meta);

fprintf('Current time: %s\n', char(datetime('now','Format','HH:mm:ss')));
mefi_exp_colec = collect_mefi_data(parameters,ic_lag,base_set,...
                                   "Modes",7,...
                                   "Max_Simulation_Time",Inf,...
                                   "Save_Data",id);

fprintf("Finished Collecting MEFI Data\n");

%% Extract MEFI Data

mefi_transtime = mefi_exp_colec.mefi_data{1}.mefi_transtime;
mefi_endtime = mefi_exp_colec.mefi_data{1}.mefi_endtime;
mefi_noc = mefi_exp_colec.mefi_data{1}.mefi_clusters;

fprintf("Finished Extracting MEFI Data\n")

%% NODE EXPERIMENT

fprintf("\nCollecting NODE Data\n");

% Experiment ID
id = "node_data_26_11_16";

% Set up record of metadata
meta = fopen(id+'_metadata.txt', 'w');
fprintf(meta,"EXPERIMENT METADATA for ID: "+id+"\n"+...
             "===================\n\n");
fprintf(meta,"Date: %s\n\n",datetime('now'));
fprintf(meta,"Summary: Simulation results for NODE across a parameter "+...
             "range localized about [0,0.1,0.55,0.7], where we "+...
             "observe a transition from a 4-cluster state to a "+...
             "2-cluster state.\n\n");
fprintf(meta,"System Parameters: Simulations were run until "+...
             "convergence (no maximum walltime), tracing 7 cluster "+...
             "modes. Initial conditions were sampled from the "+...
             "density (after normalization) function "+...
             "f(x)=1.5-cos(8pi x). Parameter range was taken from "+...
             "admissible combinations of \n"+...
             "s1 = 0.00 +                 [0]\n" + ...
             "s2 = 0.10 +        [-0.05,   0, 0.05]\n"+...
             "r1 = 0.55 + [-0.10, -0.05,   0, 0.05, 0.10]\n"+...
             "r2 = 0.70 + [-0.10, -0.05,   0, 0.05, 0.10]\n"+...
             "and the particle-number range was taken from 1e3 to 1e4 "+...
             "in increments of 1e3.");

fclose(meta);

N = 1000*(1:20);
trials = 19;

fprintf('Current time: %s\n', char(datetime('now','Format','HH:mm:ss')));
node_exp_colec = collect_node_data(parameters,generate_sample,param_space,...
                                   "Max_Simulation_Time",Inf,...
                                   "Modes",7,...
                                   "N",N,...
                                   "Save_Data",id,...
                                   "Track_Times",true,...
                                   "Trials",trials);

fprintf("Finished Collecting NODE Data\n");

%%

param_space

%%



analyze_node_data([mefi_transtime,mefi_endtime,mf_noc],exp_colec);

disp('done')