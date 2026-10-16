function [return_time, return_mesh, return_mass, return_clock, return_conv, ...
          return_end] = ....
         cont_cons_lag(ic_data, params, options)
%CONT_CONS_LAG is a conservative Lagrangian scheme for the mean-field model
%of yeast metabolism. It is built from the proof-scheme for NODE but uses
%the convergence criterion from the the proof-scheme for MF.
%
%last updated 10/16/26 by Adam Petrucci
arguments (Input)
    ic_data            % initial mesh and mass
    params             % parameters for simulation
end
arguments (Input)
    options.EndTime = Inf    % maximum simulation time
                             % Default: run until achieving convergence
    options.Update = true    % whether or not to regularly print updates
                             % Default: Deliver updates
    options.Collect = false  % data to collect
                             % Default: Returns only the first reference
                             % state and the final state
    options.Tolerance = 1e-8 % convergence tolerance
                             % Default: some preliminary runs suggest this
                             % to be a reasonable tolerance for catching
                             % metastable states
    options.Track = false    % collect and return convergence data
    options.Storage_GB = 1;  % maximum data storage
end
arguments (Output)
    return_time (1,:)   % discretized time axis of simulation
    return_mesh         % evolving mesh
    return_mass         % shift in mass corresponding to mesh
    return_clock        % total real-time for simulation
    return_conv         % convergence data
    return_end          % end by convergence or walltime
end

    % Begin timer
    timer = tic;

    %****************************
    % Process Inputs
    %****************************

    % extract initial conditions
    ic = squeeze(ic_data(1,:));
    mass_curr = squeeze(ic_data(2,:));
    num_par = length(ic);

    % collect numerical parameters
    t_final = options.EndTime;
    ud = options.Update;
    tol = options.Tolerance;
    trak = options.Track;

    %****************************
    % Change coordinate System
    %****************************
    % Change coordinate system so that r2_tilde = 1 = 0.
    % WARNING: there is approximation here

    r1_tilde = mod(params.r1 - params.r2,1);
    s1_tilde = mod(params.s1 - params.r2,1);
    s2_tilde = mod(params.s2 - params.r2,1);
    ic = mod(ic-params.r2,1);
    j = find(diff(ic)<0,1);
    mass_shift = 0;
    if ~isempty(j)
        ic = [ic(j+1:end), ic(1:j)];
        mass_curr = [mass_curr(j+1:end),mass_curr(1:j)];
        mass_shift = j;
    end

    %****************************
    % Set up Scheme
    %****************************
    % Set up initial state and instantiate counters and storage objects.

    % Instantiate storage objects according to whether user calls for data
    % at every revolution or only the final state.
    collect = options.Collect;
    if collect
        max_sav = 2*floor((floor(options.Storage_GB*2^30 / (num_par*8))-1)/2)+1;
        time = zeros([max_sav,1]);               % time vector
        mesh_data = zeros([max_sav,num_par]);    % lag state vector
        mass_data = zeros([max_sav,1]);          % mass vector
    else
        time = zeros([2,1]);                     % time vector
        mesh_data = zeros([2,num_par]);     % lag state vector
        mass_data = zeros([2,1]);                % mass vector
    end
    conv_data = zeros(size(time));               % convergence vector

    rev_c = 0;               % (completed) revolution counter
    sav_c = 1;               % save counter
    sav_f = 1;               % save frequency
    time(sav_c) = 0;         % save first time
    mass_data(sav_c) = 0;    % save first mass-shift

    full_mesh = ic;    % initial state
    tt = 0;    % initial time

    % Store state at previous revolution (for convergence criterion)
    [coll_mesh,rho,atom,mass_cell] = split_sing_cont(full_mesh,mass_curr);
    prior_mesh = coll_mesh;
    prior_rho = rho;
    prior_atom = atom;

    % Automatically use optimal timestep
    fxd_dt = s1_tilde;

    % Store initial state
    [pos_save,undo_shift] = lagrange_unshift(full_mesh,params.r2);
    mesh_data(sav_c,:) = pos_save;
    mass_data(sav_c,:) = mod(mass_shift + undo_shift,num_par);

    %****************************
    % Iterate
    %****************************

    % Give progress update (if desired).
    if ud
        fprintf("Began Simulation of Mean-Field with " + ...
                "Conservative Lagrangian scheme\n");
    end

    % initialize timestep and phantom 'particle'
    dt = fxd_dt;
    phantom = 0;

    new_rev = false;    % bool for completing a revolution
    converged = false;  % bool for achieving convergence

    % Run scheme until either 1) convergence is achieved or 2) reach
    % maximum designated wall-time.
    while ~converged && (tt < t_final)

        % If just finished a new revolution return to optimal timestep
        if new_rev
            new_rev = false;
            dt = fxd_dt;
        end

        % Split into continuous and singular components
        x_ext         = [coll_mesh(end)-1,coll_mesh,coll_mesh(1)+1];
        rho_ext       = [rho(end),rho];
        cont_mass_ext = [mass_cell(end),mass_cell];

        % Determine initial population of S
        in_ind  = find(x_ext < s1_tilde,1,'last');  % last index pre-S
        out_ind = find(x_ext < s2_tilde,1,'last');  % last index in S

        % Compute initial S population
        if in_ind == out_ind
            N0 = rho_ext(in_ind) * (s2_tilde-s1_tilde);
        else
            N0 = rho_ext(in_ind)*(x_ext(in_ind+1)-s1_tilde) + ...
                      sum(cont_mass_ext(in_ind+1:out_ind-1)) + ...
                      rho_ext(out_ind)*(s2_tilde-x_ext(out_ind));
        end
        N0 = N0 + sum(atom(in_ind : out_ind-1));

        % Construct entry/exit events relative to S
        enter_events = s1_tilde-coll_mesh;
        enter_inds = enter_events > 0 & enter_events < dt;
        enter_times = flip(enter_events(enter_inds));
        exit_events = s2_tilde-coll_mesh;
        exit_inds = exit_events > 0 & exit_events < dt;
        exit_times = flip(exit_events(exit_inds));
        
        % Compute change population at each event
        flux_changes = [rho(end),rho(1:end-1)]-rho;
        enter_flux_jumps = flip(flux_changes(enter_inds));
        exit_flux_jumps = flip(-flux_changes(exit_inds));
        enter_atom_jumps = flip(atom(enter_inds));
        exit_atom_jumps = flip(-atom(exit_inds));

        % Combing changes in population into a single, ordered vector
        [event_times, jumps] = merge_sorted(enter_times,exit_times,...
                    'AuxFun',build_aux(enter_flux_jumps,exit_flux_jumps,...
                                       enter_atom_jumps,exit_atom_jumps),...
                                       'AuxD',2);
        cont_jumps = jumps(1,:);
        atom_jumps = jumps(2,:);

        % Compute population of S
        flux0 = rho_ext(in_ind)-rho_ext(out_ind);

        flux = flux0 + [0, cumsum(cont_jumps)];
        times = [0, event_times, dt];
        Ns = N0 + [0,cumsum(flux .* diff(times) + [atom_jumps,0])];

        % Compute speed in R after each event
        speeds = 1 - params.alph*Ns;
        speeds_minus = speeds(1:end-1) - params.alph*flux.*diff(times);

        % Add boundary times and compute net displacement for particle in R
        H = [0, cumsum(0.5 *(speeds(1:end-1) + speeds_minus) ...
                          .* diff(times))];

        % Evolve the phantom particle according to three cases.
        % First case: for new revolution, dt is chosen specifially to place
        % phantom at starting position.
        % Second case: phantom is in regime of unit speed
        % Third case: phantom particle interacts with R
        if phantom <= r1_tilde - dt

            phantom = phantom + dt;

        else

            % Determine interaction of phantom with R
            phantom_entry_time = max(0,r1_tilde - phantom);
            phantom_base = max(r1_tilde, phantom);
            phantom_entry_index = find(times <= phantom_entry_time,1,'last');
            phantom_entry_diff = phantom_entry_time - times(phantom_entry_index);
            phantom_Hcorr = H(phantom_entry_index) + ...
                            speeds(phantom_entry_index) .* phantom_entry_diff - ...
                            0.5*params.alph*flux(phantom_entry_index) .* phantom_entry_diff.^2;
            phantom_target = 1 - phantom_base + phantom_Hcorr;

            % Check if phantom leaves R
            if phantom_target > H(end)  % does not leave R

                phantom = phantom_base + H(end) - phantom_Hcorr;

            else                        % does leave R

                % Since r_2 = 1, this implies a new revolution
                new_rev = true;

                % Compute exit time of phantom particle
                phantom_exit_index = find(H <= phantom_target,1,'last');
                phantom_exit_diff = phantom_target - H(phantom_exit_index);
                phantom_exit_time = times(phantom_exit_index) + ...
                                    2*phantom_exit_diff/(speeds(phantom_exit_index) + ...
                                    sqrt(speeds(phantom_exit_index)^2 - ...
                                    2*params.alph*phantom_exit_diff*flux(phantom_exit_index)));
                dt = phantom_exit_time;
                phantom = 0;

                % To accomodate the shortened timestep, the functional
                % vectors are truncated and modified.
                % The if/else here is just to handle some funky counting.
                if dt > times(phantom_exit_index)

                    tau = dt - times(phantom_exit_index);
                    times = [times(1:phantom_exit_index), dt];
                    speeds = [speeds(1:phantom_exit_index), ...
                            speeds(phantom_exit_index) - params.alph*flux(phantom_exit_index)*tau];
                    H = [H(1:phantom_exit_index), ...
                         H(phantom_exit_index) + speeds(phantom_exit_index)*tau ...
                            - 0.5*params.alph*flux(phantom_exit_index)*tau^2];
                    flux = flux(1:phantom_exit_index);

                else

                    times  = times(1:phantom_exit_index);
                    speeds = speeds(1:phantom_exit_index);
                    H      = H(1:phantom_exit_index);
                    flux   = flux(1:phantom_exit_index-1);

                end
            end
        end

        % Find particles that interact with R during step
        % In new coords, r2=1 so x<r2 is trivial
        r_array = full_mesh > r1_tilde-dt;
        r_set = full_mesh(r_array);
        r_ind = find(r_array);

        % Compute particles' starting points relative to R
        % i.e. whether they start in R or reach R during the step
        base = max(r_set, r1_tilde);

        % Compute entry times in R
        entry_times_r = zeros(size(r_set));
        enter_r = r_set < r1_tilde;
        entry_times_r(enter_r) = r1_tilde - r_set(enter_r);

        % Compute H at entry time
        % This leverages monotonicity of the state vector for speed
        Hcorr = zeros([1,length(r_set)]);
        k = 1;
        for i = flip(1:length(entry_times_r))
            while k < length(times) && times(k+1) <= entry_times_r(i)
                k = k + 1;
            end
            tau = entry_times_r(i) - times(k);
            Hcorr(i) = H(k) + speeds(k)*tau - 0.5*params.alph*flux(k)*tau^2;
        end

        % Compute displacement of particles interacting with R
        % This leverages monotonicity of the state vector for speed
        targets = 1 - base + Hcorr;
        k = 1;
        for i = flip(1:length(r_set))
            while k < length(H) && H(k+1) <= targets(i)
                k = k + 1;
            end
            if k == length(H) % particle doesn't leave R in step
                full_mesh(r_ind(i)) = base(i) + H(end) - Hcorr(i);
            else              % particle left R during step
                C = targets(i) - H(k);
                discr = speeds(k)^2 - 2*params.alph*flux(k)*C;
                exit_time = times(k) + 2*C / (speeds(k) + sqrt(discr));
                full_mesh(r_ind(i)) = dt - exit_time;
            end
        end

        % Compute displacement of all particles that do not interact with R
        full_mesh(~r_array) = full_mesh(~r_array) + dt;

        % Correct for periodicity
        full_mesh = mod(full_mesh,1);

        % Re-sort state vector to ensure monotonicity
        j = find(diff(full_mesh)<0,1);
        if ~isempty(j)
            full_mesh = [full_mesh(j+1:end), full_mesh(1:j)];
            mass_curr = [mass_curr(j+1:end),mass_curr(1:j)];
            mass_shift = mod(mass_shift + j,num_par);
        end

        [coll_mesh,rho,atom,mass_cell] = split_sing_cont(full_mesh,mass_curr);

        % Update current time
        tt = tt + dt;

        % If the current step completed a revolution, save data and test
        % for convergence.
        if new_rev

            rev_c = rev_c+1;

            % Compute difference between points of Poincare map
            change = metric_wasserstein1_lag(coll_mesh,mass_curr,...
                                             prior_mesh,mass_curr,...
                                             'Rho1',rho,'Atom1',atom,...
                                             'Rho2',prior_rho,'Atom2',prior_atom);
            prior_mesh = coll_mesh;
            prior_rho = rho;
            prior_atom = atom;

            % Check for convergence
            converged = change < tol;

            % Give progress update (if desired).
            if ud && (mod(rev_c,10) == 0)
                fprintf('Reached revolution %d at in-game time %.2f and real time %.2f\n',rev_c,tt,toc(timer));
            end

            % Update data
            if collect && (converged || (mod(rev_c,sav_f) == 0))

                sav_c = sav_c + 1;
                time(sav_c) = tt;
                [pos_save,undo_shift] = lagrange_unshift(full_mesh,params.r2);
                mesh_data(sav_c,:) = pos_save;
                mass_data(sav_c) = mod(mass_shift + undo_shift,num_par);
                if trak
                    conv_data(sav_c) = change;
                end

                % Update storage objects
                if sav_c == max_sav
                    sav_c = (max_sav+1)/2;
                    for k = 2:sav_c
                        src = 2*k-1;
                        time(k) = time(src);
                        mass_data(k) = mass_data(src);
                        mesh_data(k,:) = mesh_data(src,:);
                        conv_data(k,:) = conv_data(src,:);
                    end
                    sav_f = 2*sav_f;
                end
                
            end

        end % of revolution block

    end % of while loop

    % If only request final state, store it now
    if ~collect

        % store ending time and data according to original positions
        time(2) = tt;
        [pos_save,undo_shift] = lagrange_unshift(full_mesh,params.r2);
        mesh_data(2,:) = pos_save;
        mass_data(2) = mod(mass_shift + undo_shift,num_par);
    end

    % Stop timer
    end_time = toc(timer);

    % Return data
    return_time = time(1:sav_c);
    return_conv = conv_data(1:sav_c);
    return_mesh = mesh_data(1:sav_c,:);
    return_mass = mass_data(1:sav_c,:);
    return_clock = end_time;

    % Give progress update (if desired).
    if ud
        fprintf("Completed Simulation of Mean-Field with " + ...
                "Conservative Lagrangian scheme " + ...
                "in %f seconds ",end_time);
        if tt >= t_final
            fprintf('hitting end time wall\n');
        else
            fprintf('with convergence\n');
        end
    end

    % Check if simulation ended due to convergence or walltime
    if tt >= t_final
        return_end = 0;
    else
        return_end = 1;
    end

end % of main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION(S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,mass_out] = lagrange_unshift(d,r2)
%LAGRANGE_UNSHIFT un-does the initial shift in coordinates. It's only given
%its own function to clean up the code.
%
%last updated 10/05/26 by Adam Petrucci

    x = mod(d+r2,1);
    mass_out = 0;
    j = find(diff(x)<0,1);
    if ~isempty(j)
        x = [x(j+1:end),x(1:j)];
        mass_out = j;
    end

end

function return_data = build_aux(vec1,vec2,vec3,vec4)
%BUILD_AUX sets up the auxiliary function for merging sorted vectors. The
%function changes on every iteration according enter_jumps and exit_jumps,
%this helper produces that structure
%
%last updated 09/23/26 by Adam Petrucci

    return_data =  @(ind1,ind2,hit1,hit2) ...
                    aux_helper(ind1,ind2,hit1,hit2,vec1,vec2,vec3,vec4);
end

function return_data = aux_helper(ind1,ind2,hit1,hit2,vec1,vec2,vec3,vec4)
%AUX_HELPER is the actual auxiliary function, with more variables than the
%merge_sorted method can actually take (hence the combination with
%build_aux)
%
%last updated 10/08/26 by Adam Petrucci

    j = [0;0];

    if hit1    % response to action of first vector
        j = j + [vec1(ind1);vec3(ind1)];
    end

    if hit2    % response to action of second vector
        j = j + [vec2(ind2);vec4(ind2)];

    end

    return_data = j;

end