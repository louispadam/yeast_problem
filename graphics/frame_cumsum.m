function return_data = frame_cumsum(params,axis,options)
%FRAME_CUMSUM presents simulation data as a cumulative distribution
%
%last updated 10/05/26
arguments (Input)
    params struct           % parameters used for simulation
    axis                    % axis to format
end
arguments (Input)
    options.Regions logical = false         % show regions?
    options.Region_labels logical = false   % add regions to legend?
    options.Data = {}                       % data to plot
                                                % left boundaries
                                                % masses
    options.Meta = {}                       % meta data: struct with
                                                % name (string)
                                                % color (rgb)
    options.Title = ""                      % title of axis
    options.Legend logical = false          % include legend?
end
arguments (Output)
    return_data    % axis into which the data has been plotted
end

    %****************************
    % Collect Inputs
    %****************************

    % Required Inputs
    ax = axis;

    % Optional Inputs
    data = options.Data;
    meta = options.Meta;

    % Define Spatial Parameters
    s1 = params.s1;
    s2 = params.s2;
    r1 = params.r1;
    r2 = params.r2;

    %****************************
    % Construct Figures
    %****************************

    data_l = length(data);

    hold(ax,"on");

    % Collect data to plot
    for k = 1:data_l

        % Extract data. The circle is represented [0,1), so the data never
        % includes the right endpoint. It is added manually to clean up the
        % graph
        arr = data{k};
        arr_x = [arr(1,:),1];
        arr_y = [cumsum(arr(2,:)),1];

        % If the data does not include the left endpoint, it is added
        % manually to clean up the graph
        if arr_x(1) ~= 0
            arr_x = [0,arr_x];
            arr_y = [0,arr_y];
        end

        % Set name
        name = meta{k}.name;
        if name == ""
            name = sprintf('data %d',k);
        end

        % Set Color; default is gradient of greys
        color = meta{k}.color;
        if color == [-1,-1,-1]
            color = ([220,220,220] + ([105,105,105]-[220,220,220])*k/data_l(1))/255;
        end

        % Plot
        plot(ax,arr_x,arr_y,'linewidth',2,'DisplayName',name,'Color',color);

    end

    % Parameters for plot
    ax.XLim = [0 1];
    ax.YLim = [0,1];

    % If desired, show cutoff regions
    if options.Regions
        shade_regions(ax,params,"Region_Labels",options.Region_labels);
    end

    ylabel(ax,'Cumulative Mass');
    xlabel(ax,'Position');
    title(ax,options.Title,'Fontsize',18,'FontWeight', 'bold')
    if options.Legend
        legend show
    end

    return_data = ax;

end