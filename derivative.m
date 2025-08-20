function return_data = derivative(s1,s2,r1,r2,state)

    mass_active = sum((s1<state).*(state<s2))/length(state);
    return_data = 1+influence(mass_active).*(r1<state).*(state<r2);

end