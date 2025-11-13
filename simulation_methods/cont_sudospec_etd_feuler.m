function [return_time, return_data]=cont_sudospec_etd_feuler(initial,parameters)
%
%last updated 11/07/25 by Adam Petrucci
arguments
    initial (1,:)       % initial conditions
    parameters struct   % parameters for simulation
end

    % Collect Inputs
    ic = initial;
    params = parameters;

    %****************************
    % System Parameters
    %****************************
    L=2*pi;
    a0=params.s1*2*L-L;
    a1=params.s2*2*L-L;
    b0=params.r1*2*L-L;
    b1=params.r2*2*L-L;
    del=params.del;     % for defining characteristic function via tanh
    alpha=params.alph;  % term in linear influence
    ct=params.ct;
    msz=params.m_sz;
    ud=params.update;

    %******************************
    % Set up Fourier Transform
    %******************************
    N=length(ic); Nh=N/2;     % number of points
    dx=2*L/N;                % distance between points
    X=-L+(0:N-1)*dx; X=X.';  % vector of point positions
    kvec=fftshift(-Nh:Nh-1); kvec=kvec.'; % correction vector of positions to
                                      % match MatLab Fourier convention
    kx=kvec/2;               % ???
    kx = kx.';

    %****************************
    % Define characteristic functions
    %***************************
    dd2 = del/2; % so that tanh makes jump in approx length del
    
    c_ChiR = ct(b0,b1);
    ChiR = c_ChiR(X);
    ChiR = ChiR.';
    %ChiRf = fft(ChiR);

    c_ChiS = ct(a0,a1);
    ChiS = c_ChiS(X);
    ChiS = ChiS.';

    %***************************************************
    % Start Iterations
    %**************************************************
    h=params.dt;
    t_final = params.tfin;
    U=ic;
    Uf = fft(U);
    tt=0;

    steps = round(t_final/h + 1);
    sz = steps;
    keep = 1;
    if sz > msz
        sz = msz;
        keep = steps/msz;
    end
    time = zeros([1,sz]);
    data = zeros([sz,length(ic)]);
    data(1,:) = ic;
    
    k = 2;
    here = round(keep*k);
    
    pb = round(linspace(2,steps,20));
    n_pb = 1;

    if ud
        fprintf("Began Simulation\n");
    end

    for step = 2:steps

        l_tm = exp(-2*L*1i*kx*h).*Uf;
        %nl_tm = alpha*trapz(linspace(0,1,N),ChiS.*real(ifft(Uf)))*...
        %        conv(ChiRf,Uf,'same').*(1-exp(-2*L*1i*kx*h))/(N*2*L);
        nl_tm = alpha*trapz(linspace(-L,L,N),ChiS.*real(ifft(Uf)))*...
                fft(ChiR.*real(ifft(Uf))).*(1-exp(-2*L*1i*kx*h))/(2*L);
        Uf = l_tm + nl_tm;

        tt=tt+h;

        if step == here
            time(k) = tt;
            data(k,:) = real(ifft(Uf));
            k = k+1;
            here = round(keep*k);
            
        end

        if ud && step == pb(n_pb)
            fprintf('Simulation Progress: %3.0f%%\n',100*step/steps)
            n_pb = n_pb + 1;
        end

    end

    return_time = time;
    return_data = data;

    if ud
        fprintf('Completed Simulation\n')
    end

end