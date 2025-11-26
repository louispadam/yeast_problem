function [return_time, return_data] = particle_feuler(initial,parameters)
%FORWARD_EULER simulates the yeast NODE using a forward-euler algorithm.
%
%last updated 08/30/25 by Adam Petrucci
arguments (Input)
    initial (1,:)       % initial conditions
    parameters struct   % parameters for simulation
end
arguments (Output)
    return_time (1,:)   % discretized time axis of simulation
    return_data (:,:)   % simulation results: [time,data]
end
    
    %****************************
    % Collect Inputs
    %****************************

    ic = initial;
    params = parameters;

    % Define Temporal parameters
    dt = params.dt;
    t_final = params.tfin;
    msz=params.m_sz;
    ud=params.update;

    %****************************
    % Set up iteration
    %****************************

    % Define stepping for iteration
    steps = round(t_final/dt + 1);
    sz = steps;
    keep = 1;   % frequency with which to store iteration

    % If default time-vector is longer than permitted, replace with max
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
        d = mod(d + dt*derivative(params,d),1);
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

    % Return data
    return_time = time;
    return_data = data;

    % update if desired
    if ud
        fprintf('Completed Simulation\n')
    end

end