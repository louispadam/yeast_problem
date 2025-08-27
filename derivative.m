function return_data = derivative(params,state)

    s1 = params.s1;
    s2 = params.s2;
    r1 = params.r1;
    r2 = params.r2;

    mass_active = sum((s1<state).*(state<s2))/length(state);
    return_data = 1+influence(mass_active).*(r1<state).*(state<r2);

end