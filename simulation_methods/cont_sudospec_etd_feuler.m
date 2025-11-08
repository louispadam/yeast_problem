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
    c0=(params.r2+0.2*(1-params.r2))*2*L-L;
    c1=(params.r2+0.8*(1-params.r2))*2*L-L;
    del=params.del;     % for defining characteristic function via tanh
    eps=params.eps;     % diffusion coefficient
    alpha=params.alph;  % term in linear influence
    ct=params.ct;
    msz=params.m_sz;

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
    ChiRf = fft(ChiR);
    ChiRf = ChiRf(1,1:length(dx));

    c_ChiS = ct(a0,a1);
    ChiS = c_ChiS(X);
    ChiS = ChiS.';

    %c_ChiP = ct(c0,c1);
    %ChiP = c_ChiP(X);
    %ChiP = ChiP.';

    %***************************************************
    % Start Iterations
    %**************************************************
    h=params.dt;
    t_final = params.tfin;
    %Symb=1+h*(eps^2*kx.^2*8*pi^2+1i*kx*4*pi);
    %squash=kx.^2<(max(kx.^2)/4);
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
    for step = 2:steps
        %Uh=fft(U);
        %URh=fft(U.*ChiR);
        %PS=trapz(X,ChiS.*U)/2/L;   % interaction term (but normalized???)
        %Uh=(Uh+4*pi*h*1i*alpha*kx.*(URh*PS))./Symb;
        %Uh=Uh.*squash;    % ???
        %U=real(ifft(Uh));
        %U=U-(U<0).*(ChiP.*U); % restore U to positive in this region

        l_tm = exp(-2*L*1i*kx*h).*Uf;
        nl_tm = alpha*trapz(ChiS,real(ifft(Uf)))*conv(ChiRf,Uf).*(1-exp(-2*L*1i*kx*h))/(2*L);
        Uf = l_tm + nl_tm;

        tt=tt+h;

        if step == here
            time(k) = tt;
            data(k,:) = real(ifft(Uf));
            k = k+1;
            here = round(keep*k);
            fprintf('Done step %d of %d.\n',step,steps)
        end

        %time(step) = tt;
        %data(step,:) = U;

    end

    return_time = time;
    return_data = data;
end