function [return_time, return_data, return_clock, return_conv] = particle_proof(initial,params,options)
%PARTICLE_PROOF simulates the yeast NODE using a scheme inspired from our
%proof of the mean-field limit. It proceeds until convergence of a Poincare
%map defined by some reference particle.
%
%last updated 08/03/26 by Adam Petrucci
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
end

    % Begin timer
    tic

    %****************************
    % Collect Inputs
    %****************************

    ic = initial;

    t_final = options.EndTime;
    ud = options.Update;
    tol = options.Tolerance;
    trak = options.Track;

    %****************************
    % Change coordinate System
    %****************************
    % Change coordinate system so that r2_tilde = 1 = 0.

    r1_tilde = mod(params.r1 - params.r2,1);
    s1_tilde = mod(params.s1 - params.r2,1);
    s2_tilde = mod(params.s2 - params.r2,1);
    ic = mod(ic-params.r2,1);
    N = length(ic);      % number of particles

    %****************************
    % Set up Scheme
    %****************************
    % Set up initial state for proof-scheme, and instantiate counters and
    % storage objects.

    % Scheme requires particles to be ordered.
    % Saving labels preserves particle identities to be returned in same
    % ordering as was inputted.
    [d, labels] = sort(ic);   % iteration vector for data
    tt = 0;                   % current time

    % Find admissible reference particle.
    % If no particle in S then evolve particles (with unit speed because
    % there are no signaling particles) until leader reaches s1; leader
    % serves as reference. Otherwise, reference is first particle past s1.
    % ref_lab is label of reference particle and ref_pos is inital (after 
    % possible rotation to s1) position of reference particle.
    if ~any((d > s1_tilde) .* (d < s2_tilde))
        if ~any(d < s1_tilde)
            ref_lab = N;
            ref_pos = d(end);
            tt = s1_tilde - ref_pos + 1;
        else
            ref_lab = find(d < s1_tilde,'last');
            ref_pos = d(ref_lab);
            tt = s1_tilde - ref_pos;
        end
        d = mod(d + tt,1);
        j = find(diff(d)<0,1);
        if ~isempty(j)
            d = [d(j+1:end), d(1:j)];
            labels = [labels(j+1:end), labels(1:j)];
        end
    else
        ref_lab = find(d > s1_tilde,1,'first');
        ref_pos = d(ref_lab);
    end
    ref_lab = labels(ref_lab);

    % State upon previous revolution (for convergence criterion)
    prior = d;

    % Automatically use optimal timestep
    fxd_dt = s1_tilde;

    % Instantiate storage objects according to whether user calls for data
    % at every revolution or only the final state.
    collect = options.Collect;
    if collect
        time = zeros([1,1000]);            % time vector (dynamic)
        data = zeros([1000,length(ic)]);   % state vector (dynamic)
    else
        time = zeros(1);
        data = zeros([1,length(ic)]);
    end
    conv_data = zeros(size(time));

    % Assign first time and state
    kk = 1;                            % revolution counter
    time(kk) = tt;                     % save first time
    data(kk,labels) = d;               % save first state

    %****************************
    % Iterate
    %****************************

    % Give progress update (if desired).
    if ud
        fprintf("Began Simulation\n");
    end

    new_rev = false;    % bool for completing a revolution
    converged = false;  % bool for acheiving convergence

    % Run scheme until either 1) convergence is achieved or 2) reach
    % maximum designated simulation time.
    while ~converged && (tt < t_final)

        % Check if standard step would skip revolution (relative to
        % reference particle). If so, adjust timestep.
        dt = fxd_dt;
        check_ref = ref_pos - d(labels == ref_lab);
        if (check_ref > 0) && (check_ref < fxd_dt) && ...
                              ~(abs(check_ref) < 1e-10)
           dt = check_ref;
           new_rev = true;
        end

        % Find particles that interact with S during step
        s_set = d(d > s1_tilde - dt & d < s2_tilde);

        % Find particles that enter S during step
        enter_s = s1_tilde - s_set;
        enter_s = enter_s(enter_s > 0);

        % Find particles that leave S during step
        exit_s  = s2_tilde - s_set;
        exit_s  = exit_s(exit_s < dt);

        % Order event times and compute respective changes in S population
        event_times = [enter_s, exit_s];
        jumps = [ones(size(enter_s)),-ones(size(exit_s))];
        [event_times,idx] = sort(event_times);
        jumps = jumps(idx);

        % Compute S population after each event
        N0 = sum(d >= s1_tilde & d < s2_tilde);  % initial population
        Ns = N0 + [0, cumsum(jumps)];

        % Compute speed in R after each event
        speeds = 1 - params.alph*Ns/N;

        % Add boundary times and compute net displacement for particle in R
        times = [0, event_times, dt];
        H = [0, cumsum(speeds .* diff(times))];

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
            Hcorr(i) = H(k) + speeds(k)*(entry_times_r(i) - times(k)); % H at entry time
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
                exit_time = (targets(i)-H(k))/(speeds(k)) + times(k);
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

            kk = kk+1;

            % Update data storage
            if collect
                time(kk) = tt;
                data(kk,labels) = d; % store according to original positions
                
                % Check if storage vectors need to be extended
                if kk == length(time)
                    time(end + 1000) = 0;
                    data(end + 1000,:) = 0;
                    conv_data(end + 1000) = 0;
                end
            end

            % Compute difference between points of Poincare map
            change = metric_wasserstein1(d,prior);
            if trak
                conv_data(kk) = change;
            end

            % Check for convergence
            converged = change < tol;

            % Reset revolution counter
            new_rev = false;

            % Update previous revolution state
            prior = d;

            % Give progress update (if desired).
            if ud && (mod(kk,100) == 0)
                fprintf('Reached revolution %d at time %.2f\n',kk,tt);
            end

        end % of revolution block

    end % of while loop

    % If only request final state, store it now
    if ~collect
        time(2) = tt;
        data(2,labels) = d; % store according to original positions
    end

    % Stop timer
    end_time = toc;

    % Return data, truncating storage objects appropriately and reverting
    % the change of coordinates
    data = mod(data + params.r2,1);
    if collect
        return_time = time(1:kk);
        data = data(1:kk,:);
        return_data = data;
        return_conv = conv_data(2:kk);
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
            fprintf('hitting end time wall');
        else
            fprintf('with convergence\n');
        end
    end

end % of whole function