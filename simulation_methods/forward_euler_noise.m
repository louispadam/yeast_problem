function [return_time, return_data] = forward_euler_noise(initial,parameters)
%FORWARD_EULER_NOISE simulates the yeast NODE with diffusion using a
%forward-euler algorithm. I believe this is known as Euler-Maruyama
arguments (Input)
    initial (1,:)       % initial conditions
    parameters struct   % parameters for simulation
end
arguments (Output)
    return_time (1,:)   % discretized time axis of simulation
    return_data (:,:)   % simulation results: [time,data]
end
    
    % Collection inputs
    ic = initial;
    params = parameters;

    % Define temporal parameters
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
        data(step,:) = mod(data(step-1,:) + dt*derivative_noise(params,data(step-1,:)),1);
    end

    return_time = time;
    return_data = data;

end