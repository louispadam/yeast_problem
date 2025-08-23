function [return_time, return_data]=pseudospectral(ic,params)

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

    %******************************
    % Set up Fourier Transform
    %******************************
    n=10; N=2^n; Nh=N/2;     % number of points
    dx=2*L/N;                % distance between points
    X=-L+(0:N-1)*dx; X=X.';  % vector of point positions
    kvec=fftshift(-Nh:Nh-1); kvec=kvec.'; % correction vector of positions to
                                      % match MatLab Fourier convention
    kx=kvec/2;               % ???

    %****************************
    % Define characteristic functions
    %***************************
    ChiR = 0.25*(tanh((X-b0)/del)+1).*(tanh((b1-X)/del)+1);
    ChiS = 0.25*(tanh((X-a0)/del)+1).*(tanh((a1-X)/del)+1);
    ChiP = 0.25*(tanh((X-c0)/del)+1).*(tanh((c1-X)/del)+1); % region where we restore positivity

    %***************************************************
    % Start Iterations
    %**************************************************

    h=params.dt;
    t_final = params.tfin;
    Symb=1+h*(eps^2*kx.^2*8*pi^2+1i*kx*4*pi);
    squash=kx.^2<(max(kx.^2)/4);
    U=ic;
    tt=0;

    time = [tt];
    data = U;

    while (tt+h<t_final)
        Uh=fft(U);
        URh=fft(U.*ChiR);
        PS=trapz(X,ChiS.*U)/2/L;   % interaction term (but normalized???)
        Uh=(Uh+4*pi*h*1i*alpha*kx.*(URh*PS))./Symb;
        Uh=Uh.*squash;    % ???
        U=real(ifft(Uh));
        U=U-(U<0).*(ChiP.*U); % restore U to positive in this region

        tt=tt+h;

        time = cat(1,time,[tt]);
        data = cat(2,data,U);
    end

    return_time = time;
    return_data = data.';
end