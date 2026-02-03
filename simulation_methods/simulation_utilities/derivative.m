function return_data = derivative(parameters,curr_state)
%DERIVATIVE calculates the derivative of a particle assuming linear interaction.
%
%last updated 10/07/25 by Adam Petrucci
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
    eps = params.eps;
    dt = params.dt;
    alph = params.alph;
    ct = params.ct;

    % Construct Cutoffs
    ChiR = ct(r1,r2);
    ChiS = ct(s1,s2);

    % Calculate Derivative
    mass_active = sum(ChiS(state))/length(state);
    noise = sqrt(2*dt)*eps*(2*rand([1,length(state)])-1);

    return_data = 1+influence(mass_active,alph).*ChiR(state)+noise;

end