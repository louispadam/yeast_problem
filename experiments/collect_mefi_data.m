function return_data = collect_mefi_data(parameters,ic_c,parameter_range,options)
%COLLECT_MF_DATA simulates MF for various parameter regimes
%
%WARNING: NEEDS TO BE UPDATED FOR CONT_CONS_LAG
%
%last updated 11/15/26 by Adam Petrucci
arguments (Input)
    parameters      % parameter struct (mostly to be overwritten)
    ic_c            % continuous IC
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

    exp_colec = struct('mefi_ic',ic_c,...
                       'mefi_data',cell([size(pr,1),1]));
    exp_colec.mefi_data = cell([size(pr,1),1]);

    %****************************
    % Run Simulations
    %****************************

    time_record = fopen(save_data+"_time_record.txt",'w');
    fprintf(time_record,'s1 s2 r1 r2 type N trial runtime_seconds\n');

    % start timer
    exp_timer = tic;

    ic_mass = squeeze(ic_c(2,:));

    for ps = 1:size(pr,1) % ps for 'parameter set'

        parameters.s1 = pr(ps,1);
        parameters.s2 = pr(ps,2);
        parameters.r1 = pr(ps,3);
        parameters.r2 = pr(ps,4);

        % Run mean-field simulation
<<<<<<< Updated upstream:experiments/collect_mefi_data.m
        [time_mefi, mesh_mefi, mass_mefi, clock_mefi, ~, end_mefi] = ...
=======
        [time_mf, mesh_mf, mass_mf, clock_mf, ~, end_mf] = ...
>>>>>>> Stashed changes:experiments/collect_mf_data.m
            cont_cons_lag(ic_c,parameters,...
                          "Update",true,...
                          "Collect",true,...
                          "Track",false,...
                          "EndTime",max_simulation_time);
        data_mefi = [[0;mass_mefi],[ic_mass;mesh_mefi]];

        % Record real-time
        fprintf(time_record,'%.2f %.2f %.2f %.2f %s %d %d %.6f\n',...
                parameters.s1,parameters.s2,...
                parameters.r1,parameters.r2,...
                'mefi',0,1,clock_mefi);

        % Detect clusters of mean-field simulation
        clusters_mefi = track_clusters(time_mefi,data_mefi,...
                                       "Modes",modes,...
                                       "Type","Lagrangian",...
                                       "IC_Mass",ic_mass);

        % Extract NOC and earliest time to NOC
        [trans_mefi,stable_mefi] = ...
            analyze_clusters(time_mefi,clusters_mefi);

        % Store results
<<<<<<< Updated upstream:experiments/collect_mefi_data.m
        exp_colec.mefi_data{ps} = ...
            struct("parameter_set",[parameters.s1,...
                                    parameters.s2,...
                                    parameters.r1,...
                                    parameters.r2],...
                   "mefi_converged",end_mefi,...
                   "mefi_endtime",time_mefi(end),...
                   "mefi_fc",squeeze(data_mefi(end,:)),...
                   "mefi_transtime",trans_mefi,...
                   "mefi_clusters",stable_mefi);
=======
        exp_colec{ps} = struct('parameter_set',[parameters.s1,...
                                                parameters.s2,...
                                                parameters.r1,...
                                                parameters.r2],...
                               'mf_endstate',stable_mf,...
                               'mf_endtime',time_mf(end),...
                               'mf_enddata',squeeze(data_mf(end,:,:)),...
                               'mf_transtime',trans_mf,...
                               'mf_clusters',end_mf);
>>>>>>> Stashed changes:experiments/collect_mf_data.m

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