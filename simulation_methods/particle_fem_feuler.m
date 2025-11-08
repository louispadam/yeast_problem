function [return_time return_data] = particle_fem_feuler(x,initial,parameters)
%
%last updated 11/07/25 by Adam Petrucci

    tic;

    x = x(1:end-1);
    initial = initial(1:end-1);
    L = length(x);

    r1 = parameters.r1;
    r2 = parameters.r2;
    s1 = parameters.s1;
    s2 = parameters.s2;
    alpha = parameters.alph;
    t_final = parameters.tfin;
    dt = parameters.dt;
    msz=parameters.m_sz;

    % adjust region endpoints to align with grid
    if mod(r1,1/L) > 10^(-14)
        fprintf('Warning: r1 value of %f was misaligned with grid.\n',r1)
        fprintf('Substituted value: %f.\n',round(r1*L)/L)
        r1 = round(r1*L)/L;
    end

    if mod(r2,1/L) > 10^(-14)
        fprintf('Warning: r2 value of %f was misaligned with grid.\n',r2)
        fprintf('Substituted value: %f.\n',round(r2*L)/L)
        r2 = round(r2*L)/L;
    end

    if mod(s1,1/L) > 10^(-14)
        fprintf('Warning: s1 value of %f was misaligned with grid.\n',s1)
        fprintf('Substituted value: %f.\n',round(s1*L)/L)
        s1 = round(s1*L)/L;
    end

    if mod(s2,1/L) > 10^(-14)
        fprintf('Warning: s2 value of %f was misaligned with grid.\n',s2)
        fprintf('Substituted value: %f.\n',round(s2*L)/L)
        s2 = round(s2*L)/L;
    end

    % find indices corresponding to region endpoints
    r1_i = find(x==r1,1);
    r2_i = find(x==r2,1);
    s1_i = find(x==s1,1);
    s2_i = find(x==s2,1);

    steps = round(t_final/dt + 1);
    sz = steps;
    keep = 1;
    if sz > msz
        sz = msz;
        keep = steps/msz;
    end
    time = zeros([1,sz]);
    data = zeros([sz,length(initial)+1]);
    data(1,:) = [initial initial(1)];
    tt = 0;
    u = initial;

    k = 2;
    here = round(keep*k);

    col = zeros([1,L]);
    col(1) = 1/L;
    col(2) = 1/(6*L);
    t_part = toeplitz(col,col);

    uf_v = ones([1,L]);
    us_v = ones([1,L]);

    tpi = chol(t_part);

    for step = 2:steps

        u_minus = [u(end),u(1:end-1)];
        u_plus = [u(2:end),u(1)];
        uf = (u_minus + u)/2;
        us = (u + u_plus)/2;

        signal = -alpha*trapz(x(s1_i:s2_i),u(s1_i:s2_i));
        uf_v(r1_i+1:r2_i) = 1+signal;
        us_v(r1_i:r2_i-1) = 1+signal;

        x_part = uf_v.*uf - us_v.*us;
        x_part = x_part.';

        deriv = tpi \ (tpi' \ x_part);

        u = u + dt*deriv.';
        %u(u<0) = 0;

        tt = tt + dt;

        if step == here
            time(k) = tt;
            data(k,:) = [u u(1)];
            k = k+1;
            here = round(keep*k);
            fprintf('Done step %d of %d.\n',step,steps)
        end

    end

    end_time = toc;

    fprintf('Simulation run in %f seconds.\n',end_time)

    return_time = time;
    return_data = data;

end