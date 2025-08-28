function return_data = derivative(parameters,curr_state)
%DERIVATIVE calculates the derivative of a particle assuming linear interaction.
arguments (Input)
    parameters struct    % parameters for simulation
    curr_state (1,:)     % current state of the system
end

    % Collect Inputs
    params = parameters;
    state = curr_state;

    % Define Spatial Parameters
    s1 = params.s1;
    s2 = params.s2;
    r1 = params.r1;
    r2 = params.r2;
    alph = params.alph;

    % Calculate Derivative
    mass_active = sum((s1<state).*(state<s2))/length(state);
    return_data = 1+influence(mass_active,alph).*(r1<state).*(state<r2);

end