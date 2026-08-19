function [return_time, return_data, return_clock] = particle_proof_4vid(initial,params,options)
%PARTICLE_PROOF_4VID simulates the yeast NODE using a scheme inspired from
%our proof of the mean-field limit. It proceeds until some given endtime.
%The timestep is a given parameter and can be chosen to be less than the
%optimal one (delta). This is to help, as the name suggests, generate
%videos of NODE.
%
%last updated 07/30/26 by Adam Petrucci
arguments (Input)
    initial (1,:)       % initial conditions
    params struct       % parameters for simulation
end
arguments (Input)
    options.Timestep = 2     % timestep for simulation
                             % Default: use optimal timestep (delta)
    options.EndTime = 50     % maximum simulation time
    options.Update = true    % whether or not to regularly print updates
                             % Default: Deliver updates
    options.M_SZ = 2^15      % bound on size of storage object so Matlab
                             % doesn't complain
end
arguments (Output)
    return_time (1,:)   % discretized time axis of simulation
    return_data (:,:)   % simulation results: [time,data]
    return_clock        % total real-time for simulation
end

    % Begin timer
    tic

    %****************************
    % Collect Inputs
    %****************************

    ic = initial;

    % Define Temporal parameters
    dt = options.Timestep;
    t_final = options.EndTime;
    ud = options.Update;

    %****************************
    % Change coordinate System
    %****************************
    % Change coordinate system so that r2_tilde = 1 = 0.

    r1_tilde = mod(params.r1 - params.r2,1);
    s1_tilde = mod(params.s1 - params.r2,1);
    s2_tilde = mod(params.s2 - params.r2,1);

    % Shift initial conditions according to change in coordinates
    ic = mod(ic-params.r2,1);
    N = length(ic); % number of particles

    %****************************
    % Set up Scheme
    %****************************
    % Check admissibility of timestep and instantiate counters and
    % storage objects.

    % Check that given timestep is admissible for proof scheme. Generally,
    % this means dt < |delta|. In the adjusted coordinate system, s1 =
    % |delta|.
    if dt > s1_tilde
        fprintf('Given timestep too big;\nUsing default of %.4f instead.\n',s1_tilde)
        dt = s1_tilde;
    end

    % Scheme requires particles to be ordered.
    % Saving labels preserves particle identities to be returned in same
    % ordering as was inputted.
    [d, labels] = sort(ic);   % iteration vector for data
    tt = 0;        % iteration value for time

    % Define stepping for iteration
    steps = round(t_final/dt + 1);
    sz = steps;
    keep = 1;   % frequency with which to store iteration

    % If default time-vector would be longer than permitted, replace with 
    % given max
    msz = options.M_SZ;
    if sz > msz
        sz = msz;
        keep = steps/msz;
    end

    % Define time and space vectors to store
    time = zeros([1,sz]);
    time(1) = tt;
    data = zeros([sz,length(ic)]);
    data(1,:) = ic;

    kk = 2;      % counter for storing iteration
    here = round(keep*kk);

    % Prepare frequency of updates (if desired)
    pb = round(linspace(2,steps,20));
    n_pb = 1;

    %****************************
    % Iterate!
    %****************************

    % Give progress update (if desired).
    if ud
        fprintf("Began Simulation\n");
    end

    for step = 2:steps

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

        % Store result at previously calculated frequency
        if step == here
            time(kk) = tt;
            data(kk,labels) = d;
            kk = kk+1;
            here = round(keep*kk);
        end

        % display update if desired
        if ud && step == pb(n_pb)
            fprintf('Simulation Progress: %3.0f%%\n',100*step/steps)
            n_pb = n_pb + 1;
        end

    end

    % Stop timer
    end_time = toc;

    % Return data, reverting the change of coordinates
    return_time = time;
    return_data = mod(data + params.r2,1);
    return_clock = end_time;

    % Give progress update (if desired).
    if ud
        fprintf('Completed Simulation in %f seconds\n',end_time)
    end

end