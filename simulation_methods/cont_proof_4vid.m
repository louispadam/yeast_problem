function [return_time, return_data, return_clock] = ...
         cont_proof_4vid(initial, params, options)
%CONT_PROOF_4vid simulates the yeast MF using a scheme inspired from our
%proof of the mean-field limit. It is built on the method of
%characteristics and becomes something of an inverse problem. An effort
%was made to mirror structure of the analogous NODE scheme. This function
%method differs from cont_proof in that it takes sub-optimal timesteps for
%the purposes of generating videos.
%
%last updated 08/21/26 by Adam Petrucci
arguments (Input)
    initial (1,:)       % initial conditions
    params struct   % parameters for simulation
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
    return_clock
end

    %****************************
    % Collect Inputs
    %****************************

    % Initial data
    ic = initial;

    % Extract spatial parameters
    N = length(initial);
    dx = 1/N;
    x = linspace(0,1-dx,N);

    % Define temporal parameters
    dt = options.Timestep;
    t_final = options.EndTime;
    ud = options.Update;

    %****************************
    % Change coordinate System
    %****************************
    % Change coordinate system so that s1_tilde = 1 = 0.

    r1_tilde = mod(params.r1 - params.s1,1);
    r2_tilde = mod(params.r2 - params.s1,1);
    s2_tilde = mod(params.s2 - params.s1,1);

    % Shift initial conditions according to change in coordinates.
    % Technically, this may induce some error because s1 may not align with
    % the grid.
    [~, shift_ind] = min(abs(x-params.s1));
    d = circshift(ic,-(shift_ind-1));

    %****************************
    % Set up Scheme
    %****************************
    % Check admissibility of timestep and instantiate counters and
    % storage objects.

    % Check that given timestep is admissible for proof scheme. This means
    % dt < |delta| and dt corresponds to an integral number of spatial
    % steps (when running at unit speed). Recall s1_tilde = 1.
    opt_ts = floor((1-r2_tilde)/dx) * dx;
    if dt > opt_ts         % Requested timestep too long
        dt = opt_ts;
        fprintf(['Given timestep too big; ', ...
                 'Using default of %.6f instead.\n'],dt);
    end
    if ~(mod(dt,dx) == 0)  % Requested timestep doesn't fit spatial grid
        dt = floor(dt/dx) * dx;
        fprintf(['Timestep adjusted for spatial alignment; ', ...
                 'Using %.6f instead.\n'],dt);
    end

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

    % Compute spatial steps / time step (at unit speed)
    x_steps = floor(dt/dx);

    % Current time
    tt = 0;

    % Define time and space vectors for storage
    time = zeros([1,sz]);
    time(1) = tt;
    data = zeros([sz,length(ic)]);
    data(1,:) = ic;

    kk = 2;      % counter for storing iteration
    here = round(keep*kk);

    % Prepare frequency of updates (if desired)
    pb = round(linspace(2,steps,20));
    n_pb = 1;

    % Begin timer
    tic

    %****************************
    % Construct Stationary Objects
    %****************************
    % To save time, some objects may be computed outside the iterative loop

    % Find indices related to S
    s_set = find(x < s2_tilde);            % [0,s2)  indices

    % Find indices related to R
    r_all = find(x < r2_tilde + dt & x > r1_tilde);  % (r1,r2+dt)  indices
    r_after = find(x < r2_tilde + dt & x > r2_tilde);% (r2,r2+dt)  indices
    r_set = setdiff(r_all,r_after);                  % (r1,r2] indices
    big_R = max(r_set);      %  last index of R
    lil_R = min(r_set);      % first index of R

    % Correction matrix for quadrature over S. Once S(t0) is known, each
    % mini-step effectively translates the S region, so we remove conc at
    % ends and add conc at beginning to maintian trapezoidal integration.
    ind = length(s_set);
    corr = [-1,0,ind-1,ind] - (0:(x_steps-1)).';
    corr = mod(corr,length(x)) + 1;

    % Compute fixed quantities. Out of context, these do not all make very
    % much sense, so Ctr-F to find them in the code
    times_S = 0:dx:dt;                               % S-events
    exit_times = dt - (x(r_after) - r2_tilde);       % time to exit R
    inds_exits = floor(exit_times/dx) + 1;           % inds of the exits
    exit_diff = exit_times - times_S(inds_exits);    % for interpolation
    base = min(x(r_all), r2_tilde);                  % positioning wrt R
    exit_interp = (exit_times - (inds_exits - 1)*dx)/dx; % for interp

    %****************************
    % Iterate!
    %****************************

    % Give progress update (if desired).
    if ud
        fprintf("Began Simulation\n");
    end

    for step = 2:steps

        % compute signal during timestep
        N0 = trapz(d(s_set))*dx;
        Ns = N0 + [0,cumsum(0.5*sum([1,1,-1,-1].*d(corr),2)).']*dx;

        % compute speeds and displacements
        speeds = 1 - params.alph*Ns;
        H = [0, cumsum(0.5*(speeds(1:end-1) + speeds(2:end)) ...
                            .* diff(times_S))];

        % compute displacement upon particle exit
        Hcorr = H(end)*ones(size(r_all));
        Hcorr(r_after - min(r_all) + 1) = ...
            H(inds_exits) + speeds(inds_exits) .* exit_diff + ...
            0.5*(speeds(inds_exits+1) - speeds(inds_exits))/dx .* ...
            exit_diff.^2;

        % set target displacement according to correction for exit time
        targets = Hcorr - (base - r1_tilde);

        % I suspect I could accomplish the iteration without having to use
        % a copy of the current density, but doing so in an optimal fashion
        % may require some investment, so I will leave it for now.
        d_new = circshift(d,x_steps);

        % Reset storage vector for particle starting-points
        z0 = zeros(size(r_all));   % 'starting' locations

        k = length(times_S); % counter for H
        j = 1;               % counter for index of z0
        Ns_enter = NaN(size(r_all)); % S-pop at time of entry
        for i = 1:length(r_all)

            % Determine time of particle entry in R
            while k > 1 && H(k-1) >= targets(i)
                k = k - 1;
            end
            if k == 1 % particle started in R in this step
                z0(i) = base(i) - Hcorr(i);
            else      % particle entered R in this step
                disp_diff = targets(i) - H(k-1);
                entry_time = times_S(k-1) + 2*disp_diff/(speeds(k-1) + ...
                             sqrt(speeds(k-1)^2 + ...
                             2*(speeds(k)-speeds(k-1))/dx*disp_diff));
                z0(i) = mod(r1_tilde - entry_time,1);
                entry_interp = (entry_time-times_S(k-1)) / dx;
                Ns_enter(i) = (1-entry_interp) * Ns(k-1) + ...
                                  entry_interp * Ns(k);
            end

            % Solve inverse problem for starting position with linear
            % interpolation.
            zi = z0(i);
            if x(j) > zi
                j = 1;
            end
            while j < length(x) & x(j+1) < zi
                j = j+1;
            end
            if j == big_R          % particle is near r2 so interp breaks
                if zi > r2_tilde
                    d_new(r_all(i)) = d(j+1);
                else
                    d_new(r_all(i)) = d(j);
                end
            elseif j == lil_R - 1  % particle is near r1 so interp breaks
                if zi > r1_tilde
                    d_new(r_all(i)) = d(j+1);
                else
                    d_new(r_all(i)) = d(j);
                end
            else                   % away from boundary, interpolate
                if j == length(x)
                    d_new(r_all(i)) = d(j) + (d(1)-d(j))*(zi-x(j))/dx;
                else
                    d_new(r_all(i)) = d(j) + (d(j+1)-d(j))*(zi-x(j))/dx;
                end
            end
        end

        % save indices of those that entered
        that_enter = r_all(~isnan(Ns_enter));
        Ns_enter = Ns_enter(~isnan(Ns_enter));

        % Apply characteristic jumps
        Ns_exit = (1-exit_interp) .* Ns(inds_exits) + ...
                      exit_interp .* Ns(inds_exits+1); % S-pop at time of exit
        d_new(r_after)    = d_new(r_after) .* ...
                             exp(-params.alph * Ns_exit); % jump upon entry
        d_new(that_enter) = d_new(that_enter) .* ...
                             exp(params.alph * Ns_enter); % jump upon exit
        
        % Update time and state
        d = d_new;
        tt = tt + dt;

        % Store result at previously calculated frequency
        if step == here
            time(step) = tt;
            data(step,:) = circshift(d,shift_ind - 1);
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

    % Return data
    return_clock = end_time;
    return_time = time;
    return_data = data;

    % Give progress update (if desired).
    if ud
        fprintf('Completed Simulation in %f seconds\n',end_time)
    end

end