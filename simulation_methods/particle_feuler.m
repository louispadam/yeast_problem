function [return_time, return_data, return_clock] = particle_feuler(initial,params,options)
%FORWARD_EULER simulates the yeast NODE using a forward-euler algorithm. It
%can handle both noise (via Euler-Maruyama) and noise-less scenarios
%
%last updated 07/30/26 by Adam Petrucci
arguments (Input)
    initial (1,:)       % initial conditions
    params struct       % parameters for simulation
end
arguments (Input)
    options.Timestep = 0.01  % timestep for simulation
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

    % Begin timer
    tic
    
    %****************************
    % Collect Inputs
    %****************************

    ic = initial;

    dt = options.Timestep;
    t_final = options.EndTime;
    ud = options.Update;

    %****************************
    % Set up iteration
    %****************************

    % Define stepping for iteration
    steps = round(t_final/dt + 1);
    sz = steps;
    keep = 1;   % frequency with which to store iteration

    % If default time-vector is longer than permitted, replace with max
    msz = options.M_SZ;
    if sz > msz
        sz = msz;
        keep = steps/msz;
    end

    % Define time and space vectors to store
    time = zeros([1,sz]);
    data = zeros([sz,length(ic)]);
    data(1,:) = ic;

    d = ic;   % iteration vector for data
    tt = 0;   % iteration value for time

    k = 2;      % counter for storing iteration
    here = round(keep*k);

    % Prepare frequency of updates (if desired)
    pb = round(linspace(2,steps,20));
    n_pb = 1;

    %****************************
    % Iterate!
    %****************************

    % updates if desired
    if ud
        fprintf("Began Simulation\n");
    end

    for step = 2:steps

        % Iterate
        d = mod(d + dt*derivative(params,d,dt),1);
        tt = tt + dt;

        % Store result at previously calculated frequency
        if step == here
            time(k) = tt;
            data(k,:) = d;
            k = k+1;
            here = round(keep*k);
        end

        % display update if desired
        if ud && step == pb(n_pb)
            fprintf('Simulation Progress: %3.0f%%\n',100*step/steps)
            n_pb = n_pb + 1;
        end

    end

    end_time = toc;

    % Return data
    return_time = time;
    return_data = data;
    return_clock = end_time;

    % update if desired
    if ud
        fprintf('Completed Simulation in %f seconds\n',end_time)
    end

end