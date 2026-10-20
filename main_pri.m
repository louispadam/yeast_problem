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
ic_s = @(z) (1-cos(2*pi*k*z))+0.5;
ic_n = ic_s(x);
s_n = 0.5 * sum([(ic_n(2:end) + ic_n(1:end-1)) .* (x(2:end) - x(1:end-1)), ...
                 (ic_n(1) + ic_n(end)) * mod(x(1) - x(end),1)]);
ic = @(z) ic_s(z)/s_n;

generate_sample = construct_sampler(ic);
ic_c = ic(x);

lag_data = euler_to_lagrange(x,ic_c);
ic_mass = lag_data(2,:);


disp('Prepared Initial Data')

%% Little experiment

modes = 7;

% Run mean-field simulation
[time_mf, data_mf, ~, ~, end_mf] = ...
                         cont_cons_lag(ic_c,parameters,...
                                       "Collect",true);

% Detect clusters of mean-field simulation
clusters_mf = track_clusters(time_mf,data_mf,"Modes",modes,"Discrete",false);

% Derive cluster data
[trans_mf,stable_mf] = analyze_clusters(time_mf,clusters_mf);

%%

trials = 15;
N = 1000*(1:20);
exp_colec = cell([size(N),1]);

for n = 1:length(N)

    num_par = N(n);

    exp_colec{n} = struct('num_par',num_par,...
                          'node_endstate',zeros([trials,1]),...
                          'node_endtime',zeros([trials,1]),...
                          'node_enddata',zeros([trials,N(n)]),...
                          'node_transtime',zeros([trials,1]),...
                          'node_clusters',zeros([trials,1]));

    for t = 1:trials

        % Sample initial particle data
        ic_p = sort(generate_sample(num_par));

        % Run NODE simulation
        [time_node, data_node, ~, ~, end_node] = ...
                            particle_proof(ic_p,parameters, ...
                                           "Update",true,...
                                           "Collect",true);

        % Detect clusters of NODE simulation
        clusters_node = track_clusters(time_node,data_node,...
                                               "Modes",modes,...
                                               "Discrete",true);

        % Derive cluster data
        [trans_node,stable_node] = analyze_clusters(time_node,clusters_node);

        exp_colec{n}.node_endstate(t) = end_node;
        exp_colec{n}.node_endtime(t) = time_node(end);
        exp_colec{n}.node_enddata(t,:) = squeeze(data_node(end,:));
        exp_colec{n}.node_transtime(t) = trans_node;
        exp_colec{n}.node_clusters(t) = stable_node;

    end % end of NODE trials

end % for N

save("26_10_08_experiment.mat","exp_colec");

disp('done NODE')



