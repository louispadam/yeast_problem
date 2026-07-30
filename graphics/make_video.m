function return_data = make_video(x_data,params,fig,options)
arguments (Input)
    x_data              % discretization in x-coordinate
    params struct       % parameters used for simulation
    fig                 % figure to work from
end
arguments (Input)
    options.Frame_rate = 1                  % frame rate
    options.Pause_rate = 1/30               % pause rate
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
arguments (Output)
    return_data      % bool for success or failure
end

    %****************************
    % Collect Inputs
    %****************************

    clf(fig);

    % Required Inputs
    ax = axes(fig);

    % Optional Inputs
    reg = options.Regions;
    data = options.Data;
    meta = options.Meta;
    time = options.Time;
    ptfac = options.Frame_rate;

    %****************************
    % Run Animation
    %****************************

    v = VideoWriter('output.mp4','MPEG-4');
    v.FrameRate = 1/options.Pause_rate;
    open(v);

    to_send = cell([1,length(data)]);
    for k = 1:length(data)
        arr = data{k};
        to_send{k} = arr(1,:);
    end
    frame(x_data,params,ax, ...
          Data = to_send, ...
          Meta = meta, ...
          Title = options.Title, ...
          Legend = options.Legend, ...
          Regions = reg, ...
          Region_labels = options.Region_labels);
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

            % Collect slices of data for current frame
            b = length(data) + 2 * reg;
            m = 0;
            for k = 1:length(data)
                arr = data{k};
                if meta{k}.discrete
                    handles(b).YData = fatten_points_polynomial(x_data,arr(ind,:),meta{k}.thickness);
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

            % update height of regions to match data
            if ax.YLim < m
                ax.YLim = [0 m*1.1];
                if reg
                    handles(2).YData = 1.1*m*ChiR(x_data);
                    handles(1).YData = 1.1*m*ChiS(x_data);
                end
            end
            writeVideo(v, getframe(fig));

        end

    end

    return_data = 1;
end