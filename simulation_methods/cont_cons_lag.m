function [return_time, return_data, return_clock, return_conv, ...
          return_end] = ....
         cont_cons_lag(initial, params, options)
%CONT_CONS_LAG is a conservative Lagrangian scheme for the mean-field model
%of yeast metabolism. It is built from the proof-scheme for NODE but uses
%the convergence criterion from the the proof-scheme for MF.
%
%last updated 10/02/26 by Adam Petrucci
arguments (Input)
    initial (1,:)       % initial conditions
    params struct       % parameters for simulation
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
    options.Eulerian = false % return Eulerian mesh
end
arguments (Output)
    return_time (1,:)   % discretized time axis of simulation
    return_data         % simulation results: [time,data]
    return_clock        % total real-time for simulation
    return_conv         % convergence data
    return_end          % end by convergence or walltime
end

    % Begin timer
    timer = tic;

    %****************************
    % Process Inputs
    %****************************

    d0 = initial;                  % collect density
    N_eul = length(d0);            % number of cells (particles)
    ic = 0:1/N_eul:1-1/N_eul;      % cell boundaries  

    % mass in each cell by linear approximation of density
    mass_curr = ([d0(2:end) + d0(1:end-1), d0(1) + d0(end)] .* ...
                 [ic(2:end) - ic(1:end-1), mod(ic(1) - ic(end),1)])/2;

    % collect numerical parameters
    t_final = options.EndTime;
    ud = options.Update;
    tol = options.Tolerance;
    trak = options.Track;
    eul = options.Eulerian;

    %****************************
    % Change coordinate System
    %****************************
    % Change coordinate system so that r2_tilde = 1 = 0.
    % WARNING: there is approximation here

    r1_tilde = mod(params.r1 - params.r2,1);
    s1_tilde = mod(params.s1 - params.r2,1);
    s2_tilde = mod(params.s2 - params.r2,1);
    [~, shift_ind] = min(abs(ic-params.r2));  % to undo shift later
    ic = mod(ic-params.r2,1);
    j = find(diff(ic)<0,1);
    if ~isempty(j)
        ic = [ic(j+1:end), ic(1:j)];
        mass_curr = [mass_curr(j+1:end),mass_curr(1:j)];
    end

    %****************************
    % Set up Scheme
    %****************************
    % Set up initial state and instantiate counters and storage objects.

    % Scheme requires 'particles' to be ordered.
    % Saving labels helps monitor mass
    d = ic;
    tt = 0;                   % current time

    % Instantiate storage objects according to whether user calls for data
    % at every revolution or only the final state.
    collect = options.Collect;
    if collect
        time = zeros([1,1000]);                % time vector (dynamic)
        if eul
            data = zeros([1000,length(ic)]);   % eul state vector (dynamic)
        else
            data = zeros([1000,2,length(ic)]); % lag state vector (dynamic)
        end
    else
        time = zeros([2,1]);                % time vector (static)
        if eul
            data = zeros([2,length(ic)]);   % eul state vector(static)
        else
            data = zeros([2,2,length(ic)]); % lag state vector (static)
        end
    end
    conv_data = zeros(size(time));        % convergence vector

    rev_c = 1;                            % revolution counter
    time(rev_c) = tt;                     % save first time (0)

    % Store state at previous revolution (for convergence criterion)
    [x,rho,atom,mass_cell] = split_sing_cont(d,mass_curr);
    prior_mesh = x;
    prior_rho = rho;
    prior_atom = atom;

    % Automatically use optimal timestep
    fxd_dt = s1_tilde;

    % Store initial state
    if eul
        data(rev_c,:) = circshift(lagrange_to_euler(d,ic,mass_curr,...
                                                    "Rho",rho,...
                                                    "Atom",atom),...
                                  shift_ind-1);
    else
        [pos_save,mass_save] = lagrange_unshift(d,mass_curr,params.r2);
        data(rev_c,1,:) = pos_save;
        data(rev_c,2,:) = mass_save;
    end

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
        x_ext    = [x(end)-1,x,x(1)+1];
        rho_ext  = [rho(end),rho];
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
        enter_events = s1_tilde-x;
        enter_inds = enter_events > 0 & enter_events < dt;
        enter_times = flip(enter_events(enter_inds));
        exit_events = s2_tilde-x;
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
        r_array = d > r1_tilde-dt;
        r_set = d(r_array);
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
                d(r_ind(i)) = base(i) + H(end) - Hcorr(i);
            else              % particle left R during step
                C = targets(i) - H(k);
                discr = speeds(k)^2 - 2*params.alph*flux(k)*C;
                exit_time = times(k) + 2*C / (speeds(k) + sqrt(discr));
                d(r_ind(i)) = dt - exit_time;
            end
        end

        % Compute displacement of all particles that do not interact with R
        d(~r_array) = d(~r_array) + dt;

        % Correct for periodicity
        d = mod(d,1);

        % Resort state vector to ensure monotonicity, and update labels
        % accordingly.
        j = find(diff(d)<0,1);
        if ~isempty(j)
            d = [d(j+1:end), d(1:j)];
            mass_curr = [mass_curr(j+1:end),mass_curr(1:j)];
            %labels = [labels(j+1:end), labels(1:j)];
        end

        [x,rho,atom,mass_cell] = split_sing_cont(d,mass_curr);

        % Update current time
        tt = tt + dt;

        % If the current step completed a revolution, save data and test
        % for convergence.
        if new_rev

            rev_c = rev_c+1;

            % Update data storage
            if collect
                time(rev_c) = tt;
                if eul
                    data(rev_c,:) = circshift(lagrange_to_euler(d,ic,mass_curr),...
                                              shift_ind-1);
                else
                    [pos_save,mass_save] = lagrange_unshift(d,mass_curr,params.r2);
                    data(rev_c,1,:) = pos_save;
                    data(rev_c,2,:) = mass_save;
                end
                
                % Check if storage vectors need to be extended
                if rev_c == length(time)
                    time(end + 1000) = 0;
                    conv_data(end + 1000) = 0;
                    if eul 
                        data(end + 1000,:) = 0;
                    else
                        data(end + 1000,:,:) = 0;
                    end
                end
            end

            % Compute difference between points of Poincare map
            change = metric_wasserstein1_lag(x,mass_curr,...
                                             prior_mesh,mass_curr,...
                                             'Rho1',rho,'Atom1',atom,...
                                             'Rho2',prior_rho,'Atom2',prior_atom);
            prior_mesh = x;
            prior_rho = rho;
            prior_atom = atom;

            if trak
                conv_data(rev_c) = change;
            end

            % Check for convergence
            converged = change < tol;

            % Give progress update (if desired).
            if ud && (mod(rev_c,10) == 0)
                fprintf('Reached revolution %d at in-game time %.2f and real time %.2f\n',rev_c,tt,toc(timer));
            end

        end % of revolution block

    end % of while loop

    % If only request final state, store it now
    if ~collect

        % store ending time and data according to original positions
        time(2) = tt;
        if eul
            data(2,:) = circshift(lagrange_to_euler(d,ic,mass_curr),...
                                  shift_ind-1);
        else
            [pos_save,mass_save] = lagrange_unshift(d,mass_curr,params.r2);
            data(2,1,:) = pos_save;
            data(2,2,:) = mass_save;
        end
    end

    % Stop timer
    end_time = toc(timer);

    % Return data, truncating storage objects appropriately and reverting
    % the change of coordinates
    if collect
        return_time = time(1:rev_c);
        return_conv = conv_data(2:rev_c);
        if eul
            return_data = data(1:rev_c,:);
        else
            return_data = data(1:rev_c,:,:);
        end
    else
        return_time = time;
        return_conv = conv_data(2:end);
        return_data = data;
    end
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

function [x,mass_out] = lagrange_unshift(d,mass_in,r2)
%LAGRANGE_UNSHIFT un-does the initial shift in coordinates. It's only given
%its own function to clean up the code.
%
%last updated 10/05/26 by Adam Petrucci

    x = mod(d+r2,1);
    j = find(diff(x)<0,1);
    if ~isempty(j)
        x = [x(j+1:end),x(1:j)];
        mass_out = [mass_in(j+1:end),mass_in(1:j)];
    else
        mass_out = mass_in;
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