function [return_time, return_data, return_clock, return_conv, ...
          return_end] = ....
         cont_cons_lag(initial, params, options)
%CONT_CONS_LAG is a conservative Lagrangian scheme for the mean-field model
%of yeast metabolism. It is built from the proof-scheme for NODE but uses
%the convergence criterion from the the proof-scheme for MF.
%
%last updated 09/08/26 by Adam Petrucci
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
end
arguments (Output)
    return_time (1,:)   % discretized time axis of simulation
    return_data (:,:)   % simulation results: [time,data]
    return_clock        % total real-time for simulation
    return_conv         % convergence data
    return_end          % end by convergence or walltime
end

    % Begin timer
    tic

    %****************************
    % Process Inputs
    %****************************

    d0 = initial;          % collect density
    N = length(d0);        % number of cells (particles)
    ic = 0:1/N:1-1/N;      % cell boundaries  

    % mass in each cell by linear approximation of density
    mass = ([d0(2:end) + d0(1:end-1), d0(1) + d0(end)] .* ...
            [ic(2:end) - ic(1:end-1), mod(ic(1) - ic(end),1)])/2;

    % collect numerical parameters
    t_final = options.EndTime;
    ud = options.Update;
    tol = options.Tolerance;
    trak = options.Track;

    %****************************
    % Change coordinate System
    %****************************
    % Change coordinate system so that r2_tilde = 1 = 0.
    % WARNING: there is an approximation here, I 

    r1_tilde = mod(params.r1 - params.r2,1);
    s1_tilde = mod(params.s1 - params.r2,1);
    s2_tilde = mod(params.s2 - params.r2,1);
    [~, shift_ind] = min(abs(ic-params.r2));  % to undo shift later
    ic = mod(ic-params.r2,1);

    %****************************
    % Set up Scheme
    %****************************
    % Set up initial state and instantiate counters and storage objects.

    % Scheme requires 'particles' to be ordered.
    % Saving labels helps monitor mass
    [d, labels] = sort(ic);   % iteration vector for data
    tt = 0;                   % current time

    % Instantiate storage objects according to whether user calls for data
    % at every revolution or only the final state.
    collect = options.Collect;
    if collect
        time = zeros([1,1000]);            % time vector (dynamic)
        data = zeros([1000,length(ic)]);   % state vector (dynamic)
    else
        time = zeros([2,1]);              % time vector (initial and final)
        data = zeros([2,length(ic)]);     % time vector (initial and final)
    end
    conv_data = zeros(size(time));        % convergence vector

    rev_c = 1;                            % revolution counter
    time(rev_c) = tt;                     % save first time (0)

    % Store state at previous revolution (for convergence criterion)
    prior_mesh = d;              % previous positions
    prior_mass = mass(labels);   % previous mass assignments

    % Automatically use optimal timestep
    fxd_dt = s1_tilde;

    % Store initial state
    data(rev_c,:) = circshift(lagrange_to_euler(d,ic,mass(labels)),...
                              shift_ind-1);

    %****************************
    % Iterate
    %****************************

    % Give progress update (if desired).
    if ud
        fprintf("Began Simulation\n");
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

        % Construct entry/exit events relative to S
        enter_events = s1_tilde - d;
        enter_inds = find(enter_events > 0 & enter_events < dt);
        enter_times = enter_events(enter_events > 0 & enter_events < dt);
        exit_events = s2_tilde - d;
        exit_inds = find(exit_events > 0 & exit_events < dt);
        exit_times = exit_events(exit_events > 0 & exit_events < dt);

        mass_curr = mass(labels);
        density = mass_curr ./ ...
                  [d(2:end)-d(1:end-1), d(1)-d(end)+1];
        d_ext = [d(end)-1, d, d(1)+1];
        mass_ext = [mass_curr(end), mass_curr];
        density_ext = [density(end), density];

        in_ind = find(d_ext <= s1_tilde,1,'last');
        out_ind = find(d_ext <= s2_tilde,1,'last');

        if in_ind == out_ind
            N0 = density_ext(in_ind) * (s2_tilde - s1_tilde);
        else
            N0 = density_ext(in_ind) * (d_ext(in_ind+1) - s1_tilde) + ...
                sum(mass_ext(in_ind+1:out_ind-1)) + ...
                density_ext(out_ind) * (s2_tilde - d_ext(out_ind));
        end

        enter_times = flip(enter_times);
        enter_inds  = flip(enter_inds);

        exit_times = flip(exit_times);
        exit_inds  = flip(exit_inds);

        enter_jumps = density(mod(enter_inds-2,N) + 1) - density(enter_inds);
        exit_jumps = density(exit_inds) - density(mod(exit_inds-2,N) + 1);

        n_enter = length(enter_times);
        n_exit  = length(exit_times);

        event_times = zeros(1,n_enter+n_exit);
        jumps       = zeros(1,n_enter+n_exit);

        ie = 1;
        ix = 1;
        k  = 0;

        % Tolerance only for recognizing numerically identical times
        %time_tol = 100*eps(max(1,dt));

        while ie <= n_enter || ix <= n_exit

            % Only enter-events remain
            if ix > n_exit
                t_event = enter_times(ie);
                dq = enter_jumps(ie);
                ie = ie + 1;

            % Only exit-events remain
            elseif ie > n_enter
                t_event = exit_times(ix);
                dq = exit_jumps(ix);
                ix = ix + 1;
    
            % Enter-event occurs first
            elseif enter_times(ie) < exit_times(ix)%-time_tol
                t_event = enter_times(ie);
                dq = enter_jumps(ie);
                ie = ie + 1;

            % Exit-event occurs first
            elseif exit_times(ix) < enter_times(ie)%-time_tol
                t_event = exit_times(ix);
                dq = exit_jumps(ix);
                ix = ix + 1;

            % Simultaneous enter and exit event
            else
                t_event = 0.5 * (enter_times(ie) + exit_times(ix));
                dq = enter_jumps(ie) + exit_jumps(ix);
                ie = ie + 1;
                ix = ix + 1;
            end

            k = k + 1;
            event_times(k) = t_event;
            jumps(k) = dq;

        end

        event_times = event_times(1:k);
        jumps       = jumps(1:k);

        in_left  = find(d_ext < s1_tilde,1,'last');
        out_left = find(d_ext < s2_tilde,1,'last');
        flux0 = density_ext(in_left) - density_ext(out_left);
        flux = flux0 + [0, cumsum(jumps)];
        times = [0, event_times, dt];
        Ns = N0 + [0,cumsum(flux .* diff(times))];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Compute speed in R after each event
        speeds = 1 - params.alph*Ns;

        % Add boundary times and compute net displacement for particle in R
        times = [0, event_times, dt];
        H = [0, cumsum(0.5 *(speeds(1:end-1) + speeds(2:end)) ...
                          .* diff(times))];
%-------------------------------------------------------------------------
        % First case: for new revolution, dt is chosen specifially to place
        % phantom at starting position.
        % Second case: phantom is in regime of unit speed
        % Third case: phantom particle interacts with R
        %if new_rev
        %    phantom = 0;
        %elseif phantom < r1_tilde - dt
        if phantom <= r1_tilde - dt
            phantom = phantom + dt;
        else
            phantom_entry_time = max(0,r1_tilde - phantom);
            phantom_base = max(r1_tilde, phantom);
            phantom_entry_index = find(times <= phantom_entry_time,1,'last');
            phantom_entry_diff = phantom_entry_time - times(phantom_entry_index);
            phantom_Hcorr = H(phantom_entry_index) + ...
                            speeds(phantom_entry_index) .* phantom_entry_diff - ...
                            0.5*params.alph*flux(phantom_entry_index) .* phantom_entry_diff.^2;
            phantom_target = 1 - phantom_base + phantom_Hcorr;
            if phantom_target > H(end)
                phantom = phantom_base + H(end) - phantom_Hcorr;
            else
                new_rev = true;

                phantom_exit_index = find(H <= phantom_target,1,'last');
                phantom_exit_diff = phantom_target - H(phantom_exit_index);
                phantom_exit_time = times(phantom_exit_index) + ...
                                    2*phantom_exit_diff/(speeds(phantom_exit_index) + ...
                                    sqrt(speeds(phantom_exit_index)^2 - ...
                                    2*phantom_exit_diff*flux(phantom_exit_index)));
                dt = phantom_exit_time;
                phantom = 0;

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
%-----------------------------------------------------------------------

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
                %exit_time = (targets(i)-H(k))/(speeds(k)) + times(k);
                C = targets(i) - H(k);
                discr = speeds(k)^2 - 2*params.alph*flux(k)*C;
                exit_time = times(k) + 2*C / (speeds(k) + sqrt(discr));
                d(r_ind(i)) = 1 + dt - exit_time;
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
            labels = [labels(j+1:end), labels(1:j)];
        end

        % Update current time
        tt = tt + dt;

        % If the current step completed a revolution, save data and test
        % for convergence.
        if new_rev

            rev_c = rev_c+1;

            % Update data storage
            if collect
                time(rev_c) = tt;
                data(rev_c,:) = circshift(lagrange_to_euler(d,ic,mass(labels)),...
                                          shift_ind-1);
                
                % Check if storage vectors need to be extended
                if rev_c == length(time)
                    time(end + 1000) = 0;
                    data(end + 1000,:) = 0;
                    conv_data(end + 1000) = 0;
                end
            end

            % Compute difference between points of Poincare map
            change = metric_wasserstein1_lag(d,mass(labels),...
                                             prior_mesh,prior_mass);
            if trak
                conv_data(rev_c) = change;
            end

            % Check for convergence
            converged = change < tol;

            % Reset revolution counter
            %new_rev = false;

            % Update previous revolution state
            prior_mesh = d;
            prior_mass = mass(labels);

            % Give progress update (if desired).
            if ud && (mod(rev_c,100) == 0)
                fprintf('Reached revolution %d at time %.2f\n',rev_c,tt);
            end

        end % of revolution block

    end % of while loop

    % If only request final state, store it now
    if ~collect
        time(2) = tt;
        mass_curr = mass(labels);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Periodic extension of source mesh
        d_ext = [d(end)-1, d, d(1)+1];
        mass_ext = [mass_curr(end), mass_curr];
        cum_mass = [0, cumsum(mass_ext)];

        % Density in each Lagrangian cell
        density_lag = mass_ext ./ diff(d_ext);

        interp_mass = zeros(size(ic_ext));
        k = 1;
        for i = 1:length(ic_ext)

            while k < length(mass_ext) && ic_ext(i) >= d_ext(k+1)
                k = k + 1;
            end

            interp_mass(i) = cum_mass(k) ...
                + density_lag(k) * (ic_ext(i) - d_ext(k));
        end

        % Mass and average density on each homogeneous cell
        density_eul = diff(interp_mass) ./ diff(ic_ext);
        data(2,:) = circshift(density_eul,shift_ind-1);
                         % store according to original positions

                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    % Stop timer
    end_time = toc;

    % Return data, truncating storage objects appropriately and reverting
    % the change of coordinates
    if collect
        return_time = time(1:rev_c);
        data = data(1:rev_c,:);
        return_data = data;
        return_conv = conv_data(2:rev_c);
    else
        return_time = time;
        return_data = data;
        return_conv = conv_data(2:end);
    end
    return_clock = end_time;

    % Give progress update (if desired).
    if ud
        fprintf('Completed Simulation in %f seconds ',end_time);
        if tt >= t_final
            fprintf('hitting end time wall\n');
        else
            fprintf('with convergence\n');
        end
    end

    if tt >= t_final
        return_end = 0;
    else
        return_end = 1;
    end

end % of main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION(S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function return_data = lagrange_to_euler(lag_coord, eul_coord, mass)

    % Periodic extension (for easy vectorization)
    d_ext = [lag_coord(end)-1, lag_coord, lag_coord(1)+1];
    ic_ext = [eul_coord,1];
    mass_ext = [mass(end), mass];
    cum_mass = [0, cumsum(mass_ext)];

    % Compute density on each cell
    density_lag = mass_ext ./ diff(d_ext);

    % linear interpolation of mass in Eulerian cells, leveraging
    % monotonicity of vectors
    interp_mass = zeros(size(ic_ext));
    k = 1;
    for i = 1:length(ic_ext)
        while k < length(mass_ext) && ic_ext(i) >= d_ext(k+1)
              k = k + 1;
        end
        interp_mass(i) = cum_mass(k) ...
                + density_lag(k) * (ic_ext(i) - d_ext(k));
    end

    % density on Eulerian mesh
    return_data = diff(interp_mass) ./ diff(ic_ext);

end