function return_data = collect_node_data(...
    parameters,generate_sample,parameter_range,options)
%COLLECT_NODE_DATA simulates NODE for various parameter regimes
%
%last updated 11/15/26 by Adam Petrucci
arguments (Input)
    parameters        % parameter struct (mostly to be overwritten)
    generate_sample   % particle sampler (for initial data)
    parameter_range   % range of parameter sets to simulate
end
arguments (Input)
    options.Modes = 7                   % number of modes to monitor
    options.Max_Simulation_Time = 1000  % maximum simulation time
    options.N = [1000, 2000]            % number of particles
    options.Trials = 10                 % NODE runs
    options.Save_Data = "N"             % title to use if saving data
    options.Track_Times = true          % record simulation durations
end
arguments (Output)
    return_data     % Results of the experiment stored in cell array of
                    % structs. The first layer includes the fields
                    % 1) parameter_set - [s1,s2,r1,r2]
                    % 2) node_data - struct with results of simulation
                    % The node_data layer includes the fields
                    % 1) num_par - number of particles (integer)
                    % 2) node_endstate - converged (boolean)
                    % 3) node_endtime - final time (number)
                    % 4) node_enddata - final state (vector)
                    % 5) node_transtime - time to final NOC (number)
                    % 6) node_clusters - final NOC (integer)
end

    %****************************
    % Collect Inputs
    %****************************

    trials = options.Trials;
    N = options.N;
    save_data = options.Save_Data;

    pr = parameter_range; % pr for 'parameter range'

    modes = options.Modes;
    max_simulation_time = options.Max_Simulation_Time;

    %****************************
    % Prepare objects for Data Collection
    %****************************

    exp_colec = cell([size(pr,1),1]);

    %****************************
    % Run Simulations
    %****************************

    time_record = fopen(save_data+"_time_record.txt",'w');
    fprintf(time_record,'s1 s2 r1 r2 type N trial runtime_seconds\n');

    % start timer
    exp_timer = tic;

    for ps = 1:size(pr,1) % ps for 'parameter set'

        parameters.s1 = pr(ps,1);
        parameters.s2 = pr(ps,2);
        parameters.r1 = pr(ps,3);
        parameters.r2 = pr(ps,4);

        % Set up broad data storage
        exp_colec{ps} = struct('parameter_set',[parameters.s1,...
                                                parameters.s2,...
                                                parameters.r1,...
                                                parameters.r2],...
                               'node_data',{cell([length(N),1])});

        % Run through samples with desired number of particles
        for n = 1:length(N)

            num_par = N(n);

            % Set up specific data storage
            exp_colec{ps}.node_data{n} = ...
                          struct('num_par',num_par,...
                                 'node_converged',zeros([trials,1]),...
                                 'node_endtime',zeros([trials,1]),...
                                 'node_ic',zeros([trials,N(n)]),...
                                 'node_fc',zeros([trials,N(n)]),...
                                 'node_transtime',zeros([trials,1]),...
                                 'node_clusters',zeros([trials,1]));

            % Run desired number of trials
            for t = 1:trials

                % Sample initial particle data
                ic_p = sort(generate_sample(num_par));

                sim_timer = tic;

                % Run NODE simulation
                [time_node, data_node, ~, ~, end_node] = ...
                                    particle_proof(ic_p,parameters, ...
                                            "Update",true,...
                                            "Collect",true,...
                                            "EndTime",max_simulation_time);

                % Record real-time
                fprintf(time_record,'%.2f %.2f %.2f %.2f %s %d %d %.6f\n',...
                parameters.s1,parameters.s2,...
                parameters.r1,parameters.r2,...
                'node',N(n),t,toc(sim_timer));

                % Detect clusters of NODE simulation
                clusters_node = track_clusters(time_node,data_node,...
                                               "Modes",modes,...
                                               "Type","Discrete");

                % Extract NOC and earliest time to NOC
                [trans_node,stable_node] = ...
                    analyze_clusters(time_node,clusters_node);

                exp_colec{ps}.node_data{n}.node_converged(t) = ...
                    end_node;
                exp_colec{ps}.node_data{n}.node_endtime(t) = ...
                    time_node(end);
                exp_colec{ps}.node_data{n}.node_ic(t,:) = ...
                    ic_p;
                exp_colec{ps}.node_data{n}.node_fc(t,:) = ...
                    squeeze(data_node(end,:));
                exp_colec{ps}.node_data{n}.node_transtime(t) = ...
                    trans_node;
                exp_colec{ps}.node_data{n}.node_clusters(t) = ...
                    stable_node;

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
                save(save_data+".mat","exp_colec");
          end

    end % of experiment for all parameter sets

    fclose(time_record);

    return_data = exp_colec;

end % of function