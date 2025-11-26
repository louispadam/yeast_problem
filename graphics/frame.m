function return_data = frame(x_data,parameters,axis,options)
%FRAME Present simulation data on a provided axis handle. Accepts a cell
% array, each element of which is a vector of particles.
%
%last updated 10/07/25
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

    % Define Spatial Parameters
    s1 = params.s1;
    s2 = params.s2;
    r1 = params.r1;
    r2 = params.r2;
    ct = params.ct;

    %****************************
    % Compile Data
    %****************************
    data_l = length(data);
    to_plot = zeros([data_l,length(x)]);

    % Collect data to plot
    for k = 1:data_l
        arr = data{k};

        % if discrete data, thicken
        if meta{k}.discrete
            arr = fatten_points_polynomial(x,arr,meta{k}.thickness);
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

    % If desired, show cutoff regions
    if reg

        % Define cutoff functions
        ChiR = ct(r1,r2);
        ChiS = ct(s1,s2);

        % Construct cutoff vectors
        ChiR = 1.1*m*ChiR(x);
        ChiS = 1.1*m*ChiS(x);

        % Plot cutoffs
        temp1 = plot(ax,x,ChiR,'-','Color',[0,0,0]/255,'linewidth',1,'DisplayName','R');
        temp2 = plot(ax,x,ChiS,'-','Color',[150,150,150]/255,'linewidth',1,'DisplayName','S');

        % If not desired, remove cutoff labels from legend
        if ~regl
            set(get(get(temp1, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
            set(get(get(temp2, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        end
    end

    % Parameters for plot
    ax.XLim = [0 1];
    ax.YLim = [0 m*1.1];
    ylabel(ax,'Density');
    xlabel(ax,'Position');
    title(ax,tit,'Fontsize',18,'FontWeight', 'bold')
    if leg
        legend show
    end

    return_data = ax;

end