function return_data = collect_mf_data(parameters,ic_c,parameter_range,options)
%COLLECT_MF_DATA simulates MF for various parameter regimes
%
%last updated 10/15/26 by Adam Petrucci
arguments (Input)
    parameters      % parameter struct (mostly to be overwritten)
    ic_c            % continuous IC (as vector)
    parameter_range
end
arguments (Input)
    options.Modes = 7                   % number of modes to monitor
    options.Max_Simulation_Time = 1000  % maximum simulation time
    options.Save_Data = "N"             % title to use if saving data
    options.Track_Times = true          % record simulation durations
end
arguments (Output)
    return_data     % results of the experiment stored in cell array of
                    % structs, each with the following fields:
                    % 1) parameter_set - [s1,s2,r1,r2]
                    % 2) mf_endstate - final NOC (integer)
                    % 3) mf_enddata - final state ([time,data] vector)
                    % 4) mf_transtime - time to final NOC (number)
                    % 5) mf_clusters - evolution of clusters (matrix)
end

    %****************************
    % Collect Inputs
    %****************************

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

        % Run mean-field simulation
        [time_mf, data_mf, clock_mf, ~, end_mf] = ...
            cont_cons_lag(ic_c,parameters,...
                          "Update",true,...
                          "Collect",true,...
                          "Track",false,...
                          "EndTime",max_simulation_time);

        % Record real-time
        fprintf(time_record,'%.2f %.2f %.2f %.2f %s %d %d %.6f\n',...
                parameters.s1,parameters.s2,...
                parameters.r1,parameters.r2,...
                'mefi',0,1,clock_mf);

        % Detect clusters of mean-field simulation
        clusters_mf = track_clusters(time_mf,data_mf,...
                                     "Modes",modes,...
                                     "Discrete",false);

        % Extract NOC and earliest time to NOC
        [trans_mf,stable_mf] = analyze_clusters(time_mf,clusters_mf);

        % Store results
        exp_colec{ps} = struct('parameter_set',[parameters.s1,...
                                                parameters.s2,...
                                                parameters.r1,...
                                                parameters.r2],...
                               'mf_endstate',stable_mf,...
                               'mf_endtime',time_mf(end),...
                               'mf_enddata',squeeze(data_mf(end,:)),...
                               'mf_transtime',trans_mf,...
                               'mf_clusters',end_mf);

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