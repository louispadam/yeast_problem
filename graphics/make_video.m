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
        if strcmp(meta{k}.type,"lagrangian")
            to_send{k} = [arr(2,2:end);arr(1,2:end)];
        else % "discrete" or "eulerian"
            to_send{k} = arr(1,:);
        end
    end
    frame(x_data,params,ax, ...
          Data = to_send, ...
          Meta = meta, ...
          Title = options.Title, ...
          Legend = false, ...
          Regions = false, ...
          Region_labels = false);
    if options.Regions
        patches = shade_regions(ax,params,"Region_Labels",options.Region_labels);
        patch_r = patches(1);
        patch_s = patches(2);
    end

    if options.Legend
        legend show
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

    m=0;

    for ind = 1:length(time)

        if mod(ind, ptfac) == 0

            % Collect slices of data for current frame
            b = length(data);
            for k = 1:length(data)
                arr = data{k};
                if strcmp(meta{k}.type,"discrete")
                    handles(b).YData = fatten_points_polynomial(x_data,arr(ind,:),meta{k}.thickness);
                elseif strcmp(meta{k}.type,"eulerian")
                    handles(b).YData = arr(ind,:);
                else % strcmp(meta{k}.type,"lagrangian")
                    lag_coord = squeeze(arr(ind+1,2:end));
                    lag_mass = squeeze(circshift(arr(1,2:end),-arr(ind+1,1)));
                    handles(b).YData = lagrange_to_euler(lag_coord,x_data,lag_mass);
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

            ax.YLim = [0 m*1.1];
            yl = ylim(ax);
            patch_r.YData = [yl(1) yl(1) yl(2) yl(2)];
            patch_s.YData = [yl(1) yl(1) yl(2) yl(2)];

            writeVideo(v, getframe(fig));

        end

    end

    return_data = 1;
end