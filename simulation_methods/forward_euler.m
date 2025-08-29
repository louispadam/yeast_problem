function [return_time, return_data] = forward_euler(initial,parameters)
%FORWARD_EULER simulates the yeast NODE using a forward-euler algorithm.
arguments (Input)
    initial (1,:)       % initial conditions
    parameters struct   % parameters for simulation
end
arguments (Output)
    return_time (1,:)   % discretized time axis of simulation
    return_data (:,:)   % simulation results: [time,data]
end
    
    % Collect Inputs
    ic = initial;
    params = parameters;

    % Define Temporal parameters
    dt = params.dt;
    t_final = params.tfin;

    % Set up iteration
    steps = round(t_final/dt + 1);
    time = zeros([1,steps]);
    data = zeros([steps,length(ic)]);
    data(1,:) = ic;

    % Iterate!
    for step = 2:steps
        time(step) = time(step-1) + dt;
        data(step,:) = mod(data(step-1,:) + dt*derivative(params,data(step-1,:)),1);
    end

    return_time = time;
    return_data = data;

end