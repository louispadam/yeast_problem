function return_data = forward_euler(ic,dt,t_final,s1,s2,r1,r2)

    steps = round(t_final/dt + 1);

    time = zeros([1,steps]);
    data = zeros([length(ic),steps]);
    data(:,1) = ic;

    for step = 2:steps
        time(step) = time(step-1) + dt;
        data(:,step) = mod(data(:,step-1) + derivative(s1,s2,r1,r2,data(step-1)),1);
    end

    return_data = [time, data];

end