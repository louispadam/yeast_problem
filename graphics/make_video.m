function return_data = make_video(x_data,parameters,fig,options)
arguments (Input)
    x_data                       % discretization in x-coordinate
    parameters struct       % parameters used for simulation
    fig                     % figure to work from
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
    options.Time = []
    options.Title = ""                      % title of axis
    options.Legend logical = false          % include legend?
end

    %****************************
    % Collect Inputs
    %****************************

    clf(fig);

    % Required Inputs
    x = x_data;
    params = parameters;
    ax = axes(fig);

    % Optional Inputs
    reg = options.Regions;
    regl = options.Region_labels;
    data = options.Data;
    meta = options.Meta;
    tit = options.Title;
    leg = options.Legend;
    time = options.Time;

    %****************************
    % Define Temporal Parameters
    %****************************
    dt=parameters.dt;       % simulation time step
    tt=0;               % current time
    ptfac=parameters.fr;    % frame rate

    %****************************
    % Run Animation
    %****************************

    v = VideoWriter('output.mp4','MPEG-4');
    v.FrameRate = 30;         % Adjust for smoothness
    open(v);

    to_send = cell([1,length(data)]);
    for k = 1:length(data)
        arr = data{k};
        to_send{k} = arr(1,:);
    end
    frame(x,params,ax, ...
          Data = to_send, ...
          Meta = meta, ...
          Title = tit, ...
          Legend = leg, ...
          Regions = reg, ...
          Region_labels = regl);
    if reg
        ChiR = params.ct(params.r1,params.r2);
        ChiS = params.ct(params.s1,params.s2);
    end

    % Set up annotation for time-keeping
    if ~isempty(time)
        p = ax.Position;
        a = annotation('textbox', ...
            [p(1)+0.1*p(3),p(2)+0.9*p(4), 0.1, 0.1], ...
            'String', sprintf('Time: %d',time(1)), ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontWeight', 'bold', ...
            'FontSize', 11, ...
            'Color','k');
    end

    writeVideo(v, getframe(fig));

    handles = findobj(ax,'Type','line');

    for ind = 1:length(time)

        if mod(ind, ptfac) == 0 % I should be able to speed this up by
                                % putting it in the for loop

            %cla(ax,'reset');

            % Collect slices of data for current frame
            b = length(data) + 2 * reg;
            m = 0;
            for k = 1:length(data)
                arr = data{k};
                if meta{k}.discrete
                    handles(b).YData = fatten_points_polynomial(x,arr(ind,:),meta{k}.thickness);
                else
                    handles(b).YData = arr(ind,:);
                end
                m_temp = max(handles(b).YData,[],"all");
                if m_temp > m
                    m = m_temp;
                end
                b = b - 1;
            end

            % Update annotation tracking time
            if ~isempty(time)
                a.String = sprintf('Time: %.2f',time(ind));
            end

            if ax.YLim < m
                ax.YLim = [0 m*1.1];
                if reg
                    handles(1).YData = 1.1*m*ChiR(x);
                    handles(2).YData = 1.1*m*ChiS(x);
                end
            end
            writeVideo(v, getframe(fig));

        end
        tt=tt+dt; % Update time
    end

    return_data = 1;
end