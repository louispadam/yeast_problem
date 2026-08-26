function return_data = big_experiment(parameters,ic_c,p_sampler,options)
%BIG_EXPERIMENT simulates MF and NODE for various parameter regimes for the
%purposes of comparison
%
%last updated 08/25/26 by Adam Petrucci
arguments (Input)
    parameters      % parameter struct (mostly to be overwritten)
    ic_c            % continuous IC (as vector)
    p_sampler       % particle IC (as sampling function)
end
arguments (Input)
    options.Modes = 7                   % number of modes to monitor
    options.Max_Simulation_Time = 1000  % maximum simulation time
    options.N = [1000, 2000]            % number of particles
    options.Save_Data = "N"             % title to use if saving data
    options.Trials                      % NODE runs / MF runs

    % For the following, L2 = s2 - s1, L3 = r1 - s2, and L4 = r2 - r1.
    % *_start is the value at which the parameter should start, and
    % *_hop is the frequency for the parameter
    % In words, the experiment will involve all values *_start:*_hop:*_max
    % where *_max is set to ensure delta >= 0.1
    options.L2_start = 0.1
    options.L2_hop = 0.1
    options.L3_start = 0.1   
    options.L3_hop = 0.1
    options.L4_start = 0.1
    options.L4_hop = 0.1
end
arguments (Output)
    return_data    % results of experiment
end

    %****************************
    % Collect Inputs
    %****************************

    save_data = options.Save_Data;

    modes = options.Modes;
    max_simulation_time = options.Max_Simulation_Time;
    trials = options.Trials;

    N = options.N;
    L2_start = options.L2_start;
    L2_hop = options.L2_hop;
    L3_start = options.L3_start;
    L3_hop = options.L3_hop;
    L4_start = options.L3_start;
    L4_hop = options.L4_hop;

    %****************************
    % Run Simulations
    %****************************

    % start timer
    exp_timer = tic;

    % Update: s1, range of s2-s1, and first storage layer
    parameters.s1 = 0.0;
    L2 = L2_start:L2_hop:0.7;
    exp_colec = cell(size(L2));

    for l2 = 1:length(L2)

        % Update: s2, range of r1-s2, and second storage layer
        parameters.s2 = L2(l2);
        L3 = L3_start:L3_hop:(0.8-parameters.s2);
        exp_colec{l2} = cell(size(L3));

        for l3 = 1:length(L3)

            % Update: s1, range of r2-r1, and third storage layer
            parameters.r1 = parameters.s2 + L3(l3);
            L4 = L4_start:L4_hop:(0.9-parameters.r1);
            exp_colec{l2}{l3} = cell(size(L4));

            for l4 = 1:length(L4)

                % Update r2, and set storage for parameters, mean-field
                % data, and NODE data
                parameters.r2 = parameters.r1 + L4(l4);
                exp_colec{l2}{l3}{l4} = cell([1,3]);
                exp_colec{l2}{l3}{l4}{1} = struct('s1',parameters.s1,...
                                                  's2',parameters.s2,...
                                                  'r1',parameters.r1,...
                                                  'r2',parameters.r2);
                exp_colec{l2}{l3}{l4}{2} = cell([1,4]);
                exp_colec{l2}{l3}{l4}{3} = cell(size(N));

                % Run mean-field simulation
                [time_mf, data_mf, ~, ~, end_mf] = ...
                        cont_proof(ic_c,parameters,...
                                        "Update",false,...
                                        "Collect",true, ...
                                        "Track",false, ...
                                        "EndTime",max_simulation_time);

                % Detect clusters of mean-field simulation
                clusters_mf = track_clusters(time_mf,data_mf,...
                                                 "Modes",modes,...
                                                 "Discrete",false);

                % Store final time and state
                exp_colec{l2}{l3}{l4}{2}{1} = time_mf(end);
                exp_colec{l2}{l3}{l4}{2}{2} = squeeze(data_mf(end,:));

                % If simulation converged, store cluster data, otherwise
                % leave flags
                if end_mf
                    [exp_colec{l2}{l3}{l4}{2}{3},exp_colec{l2}{l3}{l4}{2}{4}] = ...
                        analyze_clusters(time_mf,clusters_mf);
                else
                    exp_colec{l2}{l3}{l4}{2}{3} = -1;
                    exp_colec{l2}{l3}{l4}{2}{4} = 0;
                end

                for n = 1:length(N)

                    % Set up storage for NODE data
                    exp_colec{l2}{l3}{l4}{3}{n} = cell([1,trials]);

                    for t = 1:trials

                        % Set up storage for each NODE trial
                        exp_colec{l2}{l3}{l4}{3}{n}{t} = cell([1,4]);

                        % Sample initial particle data
                        ic_p = sort(p_sampler(N(n)));

                        % Run NODE simulation
                        [time_node, data_node, ~, ~, end_node] = ...
                            particle_proof(ic_p,parameters, ...
                                            "Update",false,...
                                            "Collect",true,...
                                            "Track",false,...
                                            "EndTime",max_simulation_time);

                        % Detect clusters of NODE simulation
                        clusters_node = ...
                            track_clusters(time_node,data_node,...
                                           "Modes",modes,...
                                           "Discrete",true);

                        % Store final time and state
                        exp_colec{l2}{l3}{l4}{3}{n}{t}{1} = time_node(end);
                        exp_colec{l2}{l3}{l4}{3}{n}{t}{2} = ...
                                                 squeeze(data_node(end,:));

                        % If simulation converged, store cluster data,
                        % otherwise leave flags
                        if end_node
                            [exp_colec{l2}{l3}{l4}{3}{n}{t}{3},...
                            exp_colec{l2}{l3}{l4}{3}{n}{t}{4}] = ...
                                analyze_clusters(time_node,clusters_node);
                        else
                            exp_colec{l2}{l3}{l4}{3}{n}{t}{3} = -1;
                            exp_colec{l2}{l3}{l4}{3}{n}{t}{4} = 0;
                        end

                    end % end of NODE trials

                end % for N

                % Deliver progress update once done with a parameter set
                fprintf("Finished running parameter set:\n" + ...
                        "[s1, s2, s3, s4] = " + ...
                        "[0.0, %.2f, %.2f, %.2f]\n" + ...
                        "with time %s\n", ...
                        parameters.s2,parameters.r1,parameters.r2,...
                        string(duration(exp_timer, 'Format', 'hh:mm:ss')));

            end % for L4

        end % for L3

    end % for L2

    % Save, if requested
    if ~strcmp(save_data,"N")
        save(save_data,'exp_colec','-v7.3')
    end

    return_data = exp_colec;

end % of function