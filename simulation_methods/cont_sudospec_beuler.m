function [return_time, return_data]=cont_sudospec_beuler(initial,parameters)
%PSEUDOSPECTRAL simulates the yeast Vlasov-McKean PDE using pseudospectral
%techniques
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

    fprintf("epsilon is %d\n",parameters.eps);

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

    c_ChiS = ct(a0,a1);
    ChiS = c_ChiS(X);
    ChiS = ChiS.';

    c_ChiP = ct(c0,c1);
    ChiP = c_ChiP(X);
    ChiP = ChiP.';

    %ChiR = 0.25*(tanh((X-b0)/dd2)+1).*(tanh((b1-X)/dd2)+1);
    %ChiS = 0.25*(tanh((X-a0)/dd2)+1).*(tanh((a1-X)/dd2)+1);
    %ChiP = 0.25*(tanh((X-c0)/dd2)+1).*(tanh((c1-X)/dd2)+1); % region where we restore positivity

    %***************************************************
    % Start Iterations
    %**************************************************
    h=params.dt;
    t_final = params.tfin;
    Symb=1+h*(eps^2*kx.^2*8*pi^2+1i*kx*4*pi);
    squash=kx.^2<(max(kx.^2)/4);
    U=ic;
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
        Uh=fft(U);
        URh=fft(U.*ChiR);
        PS=trapz(X,ChiS.*U)/2/L;   % interaction term (but normalized???)
        Uh=(Uh+4*pi*h*1i*alpha*kx.*(URh*PS))./Symb;
        Uh=Uh.*squash;    % low-pass filter
        U=real(ifft(Uh));
        U=U-(U<0).*(ChiP.*U); % restore U to positive in this region

        tt=tt+h;

        if step == here
            time(k) = tt;
            data(k,:) = U;
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