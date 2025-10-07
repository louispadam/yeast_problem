function return_data = frame(x_data,parameters,axis,options)
%FRAME Fill an axis object to present simulation data
%
%last updated 10/6/25
arguments (Input)
    x_data (1,:) double     % discretization of domain
    parameters struct       % parameters used for simulation
    axis                    % axis to format
end
arguments (Input)
    options.Regions logical = false         % show regions?
    options.Region_labels logical = false   % add regions to legend?
    options.Data = {}                       % data to plot
    options.Meta = {}                       % meta data:
                                                % name
                                                % color
                                                % thickness
    options.Title = ""                      % title of axis
end

    %****************************
    % Collect Inputs
    %****************************
    x = x_data;
    params = parameters;
    ax = axis;

    reg = options.Regions;
    regl = options.Region_labels;
    d = options.Data;
    meta = options.Meta;
    tit = options.Title;

    % Define Spatial Parameters
    a0 = params.s1;
    a1 = params.s2;
    b0 = params.r1;
    b1 = params.r2;
    del = params.del;

    %****************************
    % Compile Data
    %****************************
    
    d_s = length(d);
    y = zeros([d_s,length(x)]);

    for k = 1:d_s

        arr = d{k};

        if length(arr) ~= length(x)
            arr = fatten_points_polynomial(x,arr,meta{k}.thickness);
        end

        y(k,:) = arr;

    end

    %****************************
    % Construct Figures
    %****************************

    hold(ax,"on");

    si = size(y);

    for k = 1:si(1)
        u = y(k,:);

        name = meta{k}.name;
        if name == ""
            name = sprintf('data %d',k);
        end

        color = meta{k}.color;
        if color == [-1,-1,-1]
            color = ([220,220,220] + ([105,105,105]-[220,220,220])*k/si(1))/255;
        end

        plot(ax,x,u,'linewidth',2,'DisplayName',name,'Color',color);

    end
    
    % Determine upper bound of figure
    m = max(y,[],"all");

    %ChiR = 1.1*m*0.25*(tanh(4*pi*(x-b0)/del)+1).*(tanh(4*pi*(b1-x)/del)+1);
    %ChiS = 1.1*m*0.25*(tanh(4*pi*(x-a0)/del)+1).*(tanh(4*pi*(a1-x)/del)+1);

    %if reg
    %    temp1 = plot(ax,x,ChiR,'-','Color',[192,192,192]/255,'linewidth',2,'DisplayName','responsive');
    %    temp2 = plot(ax,x,ChiS,'-','Color',[128,128,128]/255,'linewidth',2,'DisplayName','signaling');
    %    if ~regl
    %        set(get(get(temp1, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    %        set(get(get(temp2, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
    %    end
    %end

    ax.XLim = [0 1];
    ax.YLim = [0 m*1.1];
    ylabel(ax,'Density');
    xlabel(ax,'Position');
    title(ax,tit,'Fontsize',18,'FontWeight', 'bold')
    %legend show

    return_data = ax;

end