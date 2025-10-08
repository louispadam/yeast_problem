function return_data = animate(x_data,parameters,axis,options)
%ANIMATE Animate a collection of simulations. Accepts a cell array, each
%element of which is an 2d array of data: time x particle.
%
%last updated 10/07/25 by Adam Petrucci
arguments (Input)
    x_data (1,:) double     % discretization of domain
    parameters struct       % parameters used for simulation
    axis                    % axis to format
end
arguments (Input)
    options.Regions logical = false         % show regions?
    options.Region_labels logical = false   % add regions to legend?
    options.Data = {}                       % data to plot
    options.Meta = {}                       % meta data: struct with
                                                % name (string)
                                                % discrete (boolean)
                                                % color (rgb)
                                                % thickness (double)
    options.Title = ""                      % title of axis
    options.Legend logical = false          % include legend?
end

    %****************************
    % Collect Inputs
    %****************************

    % Required Inputs
    x = x_data;
    params = parameters;
    ax = axis;

    % Optional Inputs
    reg = options.Regions;
    regl = options.Region_labels;
    data = options.Data;
    meta = options.Meta;
    tit = options.Title;
    leg = options.Legend;

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

    % Set up annotation for time-keeping
    p = ax.Position;
    a = annotation('textbox', ...
        [p(1)+0.85*p(3),p(2)+0.9*p(4), 0.1, 0.1], ...
        'String', sprintf('Time: %d',0), ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', ...
        'FontSize', 11);

    while (tt<tmax)

        % Check frame rate
        if (abs((fix(ptfac*tt)-ptfac*tt))/ptfac<dt)

            cla(ax,'reset');
            
            % Set the array index of desired time
            ind = floor(tt/dt)+1;

            to_send = cell([1,length(data)]);
            for k = 1:length(data)
                arr = data{k};
                to_send{k} = arr(ind,:);
            end

            % Display frame for given time
            frame(x,params,ax, ...
                  Data = to_send, ...
                  Meta = meta, ...
                  Title = tit, ...
                  Legend = leg, ...
                  Regions = reg, ...
                  Region_labels = regl);

            a.String = sprintf('Time: %.2f',tt);
            
            pause(pr);

        end
        tt=tt+dt; % Update time
    end

    return_data = 1;
end