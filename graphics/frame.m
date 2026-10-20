function return_data = frame(x_data,params,axis,options)
%FRAME Present simulation data on a provided axis handle. Accepts a cell
% array, each element of which is a vector or discretized continuum of
% particles.
%
%last updated 10/18/25
arguments (Input)
    x_data (1,:) double     % discretization of domain
    params struct       % parameters used for simulation
    axis                    % axis to format
end
arguments (Input)
    options.Regions logical = false         % show regions?
    options.Region_labels logical = false   % add regions to legend?
    options.Data = {}                       % data to plot
    options.Meta = {}                       % meta data: struct with
                                                % name (string)
                                                % type
                                                % color (rgb)
                                                % thickness (double)
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
    x = x_data;
    ax = axis;

    % Optional Inputs
    data = options.Data;
    meta = options.Meta;

    %****************************
    % Compile Data
    %****************************
    data_l = length(data);
    to_plot = zeros([data_l,length(x)]);

    % Collect data to plot
    for k = 1:data_l
        arr = data{k};

        % if discrete data, thicken
        if strcmp(meta{k}.type,"discrete")
            arr = fatten_points_polynomial(x,arr,meta{k}.thickness);
        end
        if strcmp(meta{k}.type,"lagrangian")
            arr = lagrange_to_euler(squeeze(arr(1,:)),x,squeeze(arr(2,:)));
        end

        to_plot(k,:) = arr;
    end

    %****************************
    % Construct Figures
    %****************************

    hold(ax,"on");
    si = size(to_plot);

    for k = 1:si(1)
        u = to_plot(k,:);

        % Set name
        name = meta{k}.name;
        if name == ""
            name = sprintf('data %d',k);
        end

        % Set Color; default is gradient of greys
        color = meta{k}.color;
        if color == [-1,-1,-1]
            color = ([220,220,220] + ([105,105,105]-[220,220,220])*k/si(1))/255;
        end

        % Plot
        plot(ax,x,u,'linewidth',2,'DisplayName',name,'Color',color);

    end
    
    % Determine upper bound of figure
    m = max(to_plot,[],"all");

    % Parameters for plot
    ax.XLim = [0 1];
    ax.YLim = [0 m*1.1];

    % If desired, show regions
    if options.Regions
        shade_regions(ax,params,"Region_Labels",options.Region_labels);
    end

    ylabel(ax,'Density');
    xlabel(ax,'Position');
    title(ax,options.Title,'Fontsize',18,'FontWeight', 'bold')
    if options.Legend
        legend show
    end

    return_data = ax;

end