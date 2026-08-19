function [return_time, return_data, return_clock] = ....
         cont_sudospec_etd_feuler(initial, params, options)
%CONT_SUDOSPEC_ETD_FEULER simulates the yeast Vlasov-McKean PDE using 
%pseudospectral methods with exponentil time differencing.
%
%last updated 07/30/26 by Adam Petrucci
arguments (Input)
    initial (1,:)       % initial conditions
    params struct   % parameters for simulation
end
arguments (Input)
    options.Timestep = 0.01  % timestep for simulation
    options.EndTime = 40     % maximum simulation time
    options.Update = true    % whether or not to regularly print updates
                             % Default: Deliver updates
    options.M_SZ = 2^15      % bound on size of storage object
                             % so Matlab doesn't complain
end
arguments (Output)
    return_time (1,:)   % discretized time axis of simulation
    return_data (:,:)   % simulation results: [time,data]
    return_clock
end

    % Begin timer
    tic

    %****************************
    % Collect Inputs
    %****************************
    ic = initial;

    %****************************
    % System Parameters
    %****************************
    L=2*pi;
    a0=params.s1*2*L-L;
    a1=params.s2*2*L-L;
    b0=params.r1*2*L-L;
    b1=params.r2*2*L-L;
    eps=params.eps;
    alpha=params.alph;  % term in linear influence
    ct=params.ct;

    ud = options.Update;

    %******************************
    % Set up Fourier Transform
    %******************************
    N=length(ic); Nh=N/2;    % number of points
    dx=2*L/N;                % distance between points
    X=-L+(0:N-1)*dx; X=X.';  % vector of point positions
    kvec=fftshift(-Nh:Nh-1); % correction vector of positions to
    kvec=kvec.';             % match MatLab Fourier convention    
    kx=kvec/2;               % ???
    kx = kx.';

    %****************************
    % Define characteristic functions
    %***************************
    
    c_ChiR = ct(b0,b1);
    ChiR = c_ChiR(X);
    ChiR = ChiR.';

    c_ChiS = ct(a0,a1);
    ChiS = c_ChiS(X);
    ChiS = ChiS.';

    %***************************************************
    % Set up iteration
    %**************************************************
    h = options.Timestep;
    t_final = options.EndTime;
    U = ic;
    Uf = fft(U);
    tt=0;

    % Define stepping for iteration
    steps = round(t_final/h + 1);
    sz = steps;
    keep = 1;   % frequency with which to store iteration

    % If default time-vector is longer than permitted, replace with max
    msz = options.M_SZ;
    if sz > msz
        sz = msz;
        keep = steps/msz;
    end

    % Define time and space vectors to store
    time = zeros([1,sz]);
    data = zeros([sz,length(ic)]);
    data(1,:) = ic;
    
    k = 2;      % counter for storing iteration
    here = round(keep*k);
    
    % Prepare frequency of updates (if desired)
    pb = round(linspace(2,steps,20));
    n_pb = 1;

    %***************************************************
    % Iterate!
    %**************************************************

    % updates if desired
    if ud
        fprintf("Began Simulation\n");
    end

    b = -2*L*1i*kx - (eps*(2*L)^2)*kx.^2;

    for step = 2:steps

        % Iterate
        l_tm = exp(b*h).*Uf;
        nl_tm = alpha*trapz(linspace(-L,L,N),ChiS.*real(ifft(Uf)))*...
                fft(ChiR.*real(ifft(Uf))).*(1-exp(-2*L*1i*kx*h))/(2*L);
        % Below is an alternative method of calculating nl_tm. It doesn't
        % seem to work as well, but may be worth looking at later
        %nl_tm = alpha*trapz(linspace(0,1,N),ChiS.*real(ifft(Uf)))*...
        %        conv(ChiRf,Uf,'same').*(1-exp(-2*L*1i*kx*h))/(N*2*L);
        
        Uf = l_tm + nl_tm;
        tt=tt+h;

        % Store result at previously calculated frequency
        if step == here
            time(k) = tt;
            data(k,:) = real(ifft(Uf));
            k = k+1;
            here = round(keep*k);
        end

        % display update if desired
        if ud && step == pb(n_pb)
            fprintf('Simulation Progress: %3.0f%%\n',100*step/steps)
            n_pb = n_pb + 1;
        end

    end
    
    end_time = toc;

    % Return data
    return_time = time;
    return_data = data;
    return_clock = end_time;

    % update if desired
    if ud
        fprintf('Completed Simulation in %f seconds\n',end_time)
    end

end
