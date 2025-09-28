function return_data = animate(x_data,parameters,fig_num,options)
%ANIMATE Animate a collection of simulations.
%
%last updated 09/24/25 by Adam Petrucci
arguments (Input)
    x_data (1,:) double     % discretization of domain
    %y_data (:,:,:) double   % results of simulation [exp_num,time,distribution]
    parameters struct       % parameters used for simulation
    fig_num                 % figure to plot in
end
arguments (Input)
    options.Names (1,:) string = []         % names of FIRST FEW curves
    options.Colors (:,3) double = []        % colors of FIRST FEW curves
    options.Regions logical = false         % show regions?
    options.Region_labels logical = false   % add regions to legend?
    %options.Continuous = []
    %options.Discrete = []
    options.Data = {}
    options.Thickness = []
end

    %****************************
    % Collect Inputs
    %****************************
    x = x_data;
    %y = y_data;
    params = parameters;
    fn = fig_num;

    names = options.Names;
    colors = options.Colors;
    reg = options.Regions;
    regl = options.Region_labels;
    %ps = options.Continuous;
    %pp = options.Discrete;
    pt = options.Thickness;
    d = options.Data;

    %****************************
    % Define Temporal Parameters
    %****************************
    dt=params.dt;       % simulation time step
    tmax = params.tfin; % ending time
    tt=0;               % current time
    ptfac=params.fr;    % frame rate
    pr = params.pr;     % pause rate

    %****************************
    % Run Animation
    %****************************
    while (tt<tmax)

        % Check frame rate
        if (abs((fix(ptfac*tt)-ptfac*tt))/ptfac<dt)
            
            % Set the array index of desired time
            ind = floor(tt/dt)+1;

            ts = cellfun(@(x) x(ind,:),d,'UniformOutput',false);

            % Display frame for given time
            frame(x,params,fn,Names=names,Colors=colors,Index=ind, ...
                  Regions=reg,Region_labels=regl,Time=tt, ...
                  Data=ts,Thickness=pt);
            pause(pr)

        end
        tt=tt+dt; % Update time
    end

    return_data = 1;
end