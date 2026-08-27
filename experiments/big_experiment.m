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
    options.Trials = 10                 % NODE runs / MF runs

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
    % Prepare Error Collection
    %****************************

    error_data_mf = struct('parameters', {}, ...
                           'message', {}, ...
                           'identifier', {}, ...
                           'stack', {}, ...
                           'report', {});

    error_data_node = struct('parameters', {}, ...
                             'message', {}, ...
                             'identifier', {}, ...
                             'stack', {}, ...
                             'report', {});

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
                exp_colec{l2}{l3}{l4}{2} = cell([1,5]);
                exp_colec{l2}{l3}{l4}{3} = cell(size(N));

                % If simulation returns error, save parameters for
                % later investigation without terminating
                % experiment
                try

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

                    % Derive cluster data
                    [trans_mf,stable_mf] = ...
                        analyze_clusters(time_mf,clusters_mf);

                catch ME

                    fprintf("ERROR simulating MF with" + ...
                            "parameter set:\n" + ...
                            "[s1, s2, r1, r2] = " + ...
                            "[0.0, %.2f, %.2f, %.2f]\n",...
                            parameters.s2,parameters.r1,parameters.r2);
                    fprintf('%s\n', ME.message);

                    num_errors = numel(error_data_mf) + 1;

                    error_data_mf(num_errors).parameters = parameters;
                    error_data_mf(num_errors).message = ME.message;
                    error_data_mf(num_errors).identifier = ME.identifier;
                    error_data_mf(num_errors).stack = ME.stack;
                    error_data_mf(num_errors).report = getReport(ME);

                    time_mf = 0;
                    data_mf = ic_c;
                    end_mf = 0;
                    trans_mf = 0;
                    stable_mf = 0;

                end

                % Store data
                exp_colec{l2}{l3}{l4}{2}{1} = end_mf;
                exp_colec{l2}{l3}{l4}{2}{2} = time_mf(end);
                exp_colec{l2}{l3}{l4}{2}{3} = squeeze(data_mf(end,:));
                exp_colec{l2}{l3}{l4}{2}{4} = trans_mf;
                exp_colec{l2}{l3}{l4}{2}{5} = stable_mf;

                for n = 1:length(N)

                    % Set up storage for NODE data
                    exp_colec{l2}{l3}{l4}{3}{n} = cell([1,trials]);

                    for t = 1:trials

                        % Set up storage for each NODE trial
                        exp_colec{l2}{l3}{l4}{3}{n}{t} = cell([1,5]);

                        % Sample initial particle data
                        ic_p = sort(p_sampler(N(n)));

                        % If simulation returns error, save parameters for
                        % later investigation without terminating
                        % experiment
                        try

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

                            % Derive cluster data
                            [trans_node,stable_node] = ...
                                analyze_clusters(time_node,clusters_node);

                        catch ME

                            fprintf("ERROR simulating NODE with" + ...
                                "parameter set:\n" + ...
                                "[s1, s2, r1, r2] = " + ...
                                "[0.0, %.2f, %.2f, %.2f]\n",...
                                parameters.s2,parameters.r1,parameters.r2);
                            fprintf('%s\n', ME.message);

                            num_errors = numel(error_data_node) + 1;

                            error_data_node(num_errors).parameters = parameters;
                            error_data_node(num_errors).message = ME.message;
                            error_data_node(num_errors).identifier = ME.identifier;
                            error_data_node(num_errors).stack = ME.stack;
                            error_data_node(num_errors).report = getReport(ME);

                            time_node = 0;
                            data_node = ic_p;
                            end_node = 0;
                            trans_node = 0;
                            stable_node = 0;

                        end

                        % Store data
                        exp_colec{l2}{l3}{l4}{3}{n}{t}{1} = end_node;
                        exp_colec{l2}{l3}{l4}{3}{n}{t}{2} = time_node(end);
                        exp_colec{l2}{l3}{l4}{3}{n}{t}{3} = ...
                                                 squeeze(data_node(end,:));
                        exp_colec{l2}{l3}{l4}{3}{n}{t}{4} = trans_node;
                        exp_colec{l2}{l3}{l4}{3}{n}{t}{5} = stable_node;

                    end % end of NODE trials

                end % for N

                % Deliver progress update once done with a parameter set
                so_far_time = toc(exp_timer);
                fprintf("Finished running parameter set:\n" + ...
                        "[s1, s2, r1, r2] = " + ...
                        "[0.0, %.2f, %.2f, %.2f]\n" + ...
                        "with time %02d:%02d:%02d\n", ...
                        parameters.s2,parameters.r1,parameters.r2,...
                        floor(so_far_time/3600),...
                        floor(mod(so_far_time,3600)/60),...
                        floor(mod(so_far_time,60)));

                % Save, if requested
                if ~strcmp(save_data,"N")
                    save(save_data,"exp_colec");
                    save("error_mf_" + save_data, "error_data_mf");
                    save("error_node_" + save_data, "error_data_node");
                end

            end % for L4

        end % for L3

    end % for L2

    return_data = exp_colec;

end % of function