function [return_time, return_data, return_clock, return_conv, ...
          return_end] = ...
         cont_proof(initial, params, options)
%CONT_PROOF simulates the yeast MF using a scheme inspired from our
%proof of the mean-field limit. It is built on the method of
%characteristics and becomes something of an inverse problem. An effort
%was made to mirror structure of the analogous
%NODE scheme.
%
%last updated 08/25/26 by Adam Petrucci
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
    options.Tolerance = 1/length(initial) % convergence tolerance
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
    t_final = options.EndTime;
    ud = options.Update;
    tol = options.Tolerance;
    trak = options.Track;

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

    % Compute optimal timestep (that respects spatial mesh)
    fxd_dt = floor((1-r2_tilde)/dx) * dx;

    % Compute spatial steps / time step (at unit speed)
    fxd_x_steps = floor(fxd_dt/dx);

    % Current time
    tt = 0;

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
    rev_c = 1;                            % revolution counter
    time(rev_c) = tt;                     % save first time
    data(rev_c,:) = d;               % save first state

    % State upon previous revolution (for convergence criterion)
    prior = d;

    % Begin timer
    tic

    %****************************
    % Construct Stationary Objects
    %****************************
    % To save time, some objects may be computed outside the iterative loop

    % Find indices related to S
    s_set = find(x < s2_tilde);            % [0,s2)  indices

    % Find indices related to R
    fxd_r_all = find(x < r2_tilde + fxd_dt & x > r1_tilde);  % (r1,r2+dt)  inds
    fxd_r_after = find(x < r2_tilde + fxd_dt & x > r2_tilde);% (r2,r2+dt)  inds
    r_set = setdiff(fxd_r_all,fxd_r_after);              % (r1,r2] inds
    big_R = max(r_set);      %  last index of R
    lil_R = min(r_set);      % first index of R

    % Correction matrix for quadrature over S. Once S(t0) is known, each
    % mini-step effectively translates the S region, so we remove conc at
    % ends and add conc at beginning to maintian trapezoidal integration.
    ind_s2 = length(s_set);
    fxd_corr = [-1,0,ind_s2-1,ind_s2] - (0:(fxd_x_steps-1)).';
    fxd_corr = mod(fxd_corr,length(x)) + 1;

    % Compute fixed quantities. Out of context, these do not all make very
    % much sense, so Ctr-F to find them in the code
    fxd_times_S = 0:dx:fxd_dt;                               % S-events
    fxd_exit_times = fxd_dt - (x(fxd_r_after) - r2_tilde);       % time to exit R
    fxd_inds_exits = floor(fxd_exit_times/dx) + 1;           % inds of the exits
    fxd_exit_diff = fxd_exit_times - fxd_times_S(fxd_inds_exits);    % for interpolation
    fxd_base = min(x(fxd_r_all), r2_tilde);                  % positioning wrt R
    fxd_exit_interp = (fxd_exit_times - (fxd_inds_exits - 1)*dx)/dx; % for interp

    %****************************
    % Iterate!
    %****************************

    % Give progress update (if desired).
    if ud
        fprintf("Began Simulation\n");
    end

    dt = fxd_dt;
    phantom = 0;

    new_rev = false;    % bool for completing a revolution
    converged = false;  % bool for acheiving convergence

    while ~converged && (tt < t_final)

        % If just finished a new revolution or started the simulation, set
        % all variables according to optimal timestep
        if new_rev || rev_c == 1

            new_rev = false;
            dt = fxd_dt;

            x_steps = fxd_x_steps;
            r_all = fxd_r_all;
            r_after = fxd_r_after;
            corr = fxd_corr;
            times_S = fxd_times_S;
            inds_exits = fxd_inds_exits;
            exit_diff = fxd_exit_diff;
            base = fxd_base;
            exit_interp = fxd_exit_interp;

        end

        % If optimal timestep would push phantom particle past its
        % reference point (indicating a new revolution) set all variables
        % according to timestep that lands revolution
        if phantom + dt > 1

            new_rev = true;
            dt = floor((1-phantom)/dx) * dx;
            x_steps = floor(dt/dx);
            
            r_all = find(x < r2_tilde + dt & x > r1_tilde);  % (r1,r2+dt)  inds
            r_after = find(x < r2_tilde + dt & x > r2_tilde);% (r2,r2+dt)  inds
            corr = [-1,0,ind_s2-1,ind_s2] - (0:(x_steps-1)).';
            corr = mod(corr,length(x)) + 1;
            times_S = 0:dx:dt;                               % S-events
            exit_times = dt - (x(r_after) - r2_tilde);       % time to exit R
            inds_exits = floor(exit_times/dx) + 1;           % inds of the exits
            exit_diff = exit_times - times_S(inds_exits);    % for interpolation
            base = min(x(r_all), r2_tilde);                  % positioning wrt R
            exit_interp = (exit_times - (inds_exits - 1)*dx)/dx; % for interp

        end

        % compute signal during timestep
        N0 = trapz(d(s_set))*dx;
        Ns = N0 + [0,cumsum(0.5*sum([1,1,-1,-1].*d(corr),2)).']*dx;

        % compute speeds and displacements
        speeds = 1 - params.alph*Ns;
        H = [0, cumsum(0.5*(speeds(1:end-1) + speeds(2:end)) ...
                            .* diff(times_S))];

        % First case: for new revolution, dt is chosen specifially to place
        % phantom at starting position.
        % Second case: phantom is in regime of unit speed
        % Third case: phantom particle interacts with R
        if new_rev
            phantom = 0;
        elseif phantom < r1_tilde - dt || phantom > r2_tilde
            phantom = phantom + dt;
        else
            % some of these don't actually have to be recalculated. there
            % is potential for optimization
            phantom_entry_time = max(0,r1_tilde - phantom);
            phantom_base = max(r1_tilde, phantom);
            phantom_entry_index = floor(phantom_entry_time/dx) + 1;
            phantom_entry_diff = phantom_entry_time - times_S(phantom_entry_index);
            phantom_Hcorr = H(phantom_entry_index) + ...
                            speeds(phantom_entry_index) .* phantom_entry_diff + ...
                            0.5*(speeds(phantom_entry_index+1) - ...
                            speeds(phantom_entry_index))/dx .* phantom_entry_diff.^2;
            phantom_target = r2_tilde - phantom_base + phantom_Hcorr;
            phantom_exit_index = find(H <= phantom_target,1,'last');
            if phantom_exit_index == length(H)
                phantom = phantom_base + H(end) - phantom_Hcorr;
            else
                phantom_exit_diff = phantom_target - H(phantom_exit_index);
                phantom_exit_time = times_S(phantom_exit_index) + ...
                                    2*phantom_exit_diff/(speeds(phantom_exit_index) + ...
                                    sqrt(speeds(phantom_exit_index)^2 + ...
                                    2*phantom_exit_diff*(speeds(phantom_exit_index+1)- ...
                                                         speeds(phantom_exit_index))/dx));
                phantom = r2_tilde + dt - phantom_exit_time;
            end
        end

        % compute displacement upon particle exit
        Hcorr = H(end)*ones(size(r_all));
        Hcorr(r_after - min(r_all) + 1) = ...
            H(inds_exits) + speeds(inds_exits) .* exit_diff ...
                       + 0.5*(speeds(inds_exits+1) - speeds(inds_exits))/dx ...
                             .* exit_diff.^2;

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
            % If r1_tilde < dt then initial points will wrap around, and j
            % must be reset
            if x(j) > zi
                j = 1;
            end
            while j < length(x) & x(j+1) < zi
                j = j + 1;
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

        % If the current step completed a revolution, save data and test
        % for convergence.
        if new_rev

            rev_c = rev_c+1;

            % Update data storage
            if collect
                time(rev_c) = tt;
                data(rev_c,:) = circshift(d,shift_ind - 1);
                
                % Check if storage vectors need to be extended
                if rev_c == length(time)
                    time(end + 1000) = 0;
                    data(end + 1000,:) = 0;
                    conv_data(end + 1000) = 0;
                end
            end

            % Compute difference between points of Poincare map
            change = metric_wasserstein1_cont(d,prior);
            if trak
                conv_data(rev_c) = change;
            end

            % Check for convergence
            converged = change < tol;

            % Update previous revolution state
            prior = d;

            % Give progress update (if desired).
            if ud && (mod(rev_c,100) == 0)
                fprintf('Reached revolution %d\n',rev_c);
            end

        end % of revolution block

    end

    % If only request final state, store it now
    if ~collect
        time(2) = tt;
        data(2,:) = circshift(d,shift_ind - 1); % store according to original positions
    end

    % Stop timer
    end_time = toc;

    % Return data, truncating storage objects appropriately and reverting
    % the change of coordinates
    if collect
        return_time = time(1:rev_c);
        return_data = data(1:rev_c,:);
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

end