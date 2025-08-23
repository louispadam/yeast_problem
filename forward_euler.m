function [return_time, return_data] = forward_euler(ic,params)

    s1 = params.s1;
    s2 = params.s2;
    r1 = params.r1;
    r2 = params.r2;

    dt = params.dt;
    t_final = params.tfin;

    steps = round(t_final/dt + 1);

    time = zeros([1,steps]);
    data = zeros([length(ic),steps]);
    data(:,1) = ic;

    for step = 2:steps
        time(step) = time(step-1) + dt;
        data(:,step) = mod(data(:,step-1) + dt*derivative(s1,s2,r1,r2,data(:,step-1)),1);
    end

    return_time = time;
    return_data = data;

end