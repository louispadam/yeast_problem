function return_data = frame(x_data,parameters,fig_num,axis,options)
%FRAME Present a collection of simulations.
%
%last updated 09/27/25
arguments (Input)
    x_data (1,:) double     % discretization of domain
    parameters struct       % parameters used for simulation
    fig_num                 % figure to plot in
    axis
end
arguments (Input)
    options.Names (1,:) string = []         % names of FIRST FEW curves
    options.Colors (:,3) double = []        % colors of FIRST FEW curves
    options.Regions logical = false         % show regions?
    options.Region_labels logical = false   % add regions to legend?
    options.Time = 0                        % time point of simulation
    options.Index = 1                       % time step of simulation
    options.Data = {}
    options.Meta = {}
    %options.Continuous = {}   % results of simulation [exp_num,time,data]
    %options.Discrete = {}     % results of simulation [exp_num,time,data]
    options.Thickness = []
    options.Title = ""
end

    %****************************
    % Collect Inputs
    %****************************
    x = x_data;
    params = parameters;
    fn = fig_num;
    ax = axis;

    %names = options.Names;
    %colors = options.Colors;
    reg = options.Regions;
    regl = options.Region_labels;
    tt = options.Time;
    %ind = options.Index;
    %ps = options.Continuous;
    %pp = options.Discrete;
    %pt = options.Thickness;
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
    %ps_s = size(ps);
    %pp_s = size(pp);

    %y = cell([pp_s+ps_s,1,length(x)]);
    y = zeros([d_s,length(x)]);

    for k = 1:d_s

        arr = d{k};

        if length(arr) ~= length(x)
            arr = fatten_points_polynomial(x,arr,meta{k}.thickness);
        end

        y(k,:) = arr;

    end

    %if ps_s > 0
    %for i = 1:cs_s
    %    y(i,:) = ps{i};
    %end

    %for k = 1:cp_s
    %for k = 1:pp_s
    %    arr = pp{k};
    %    y(cs_s+k,:) = fatten_points_polynomial(x,arr(ind,:),0.1);
    %    y(cs_s+k,1,:) = fatten_points_polynomial(x,reshape(pp(k,ind,:),[1,pp_s(3)]),pt(k));
    %end

    %****************************
    % Construct Figures
    %****************************
    %figure(fn)
    %clf
    hold(ax,"on");

    si = size(y);
    %ci = size(colors);
    for k = 1:si(1)
        u = y(k,:);
        %u = reshape(y(k,1,:),[1,si(3)]);
        %u = squeeze(y(k,ind,:));

        name = meta{k}.name;
        if name == ""
            name = sprintf('data %d',k);
        end

        color = meta{k}.color;
        if color == [-1,-1,-1]
            color = ([220,220,220] + ([105,105,105]-[220,220,220])*k/si(1))/255;
        end

        plot(ax,x,u,'linewidth',2,'DisplayName',name,'Color',color);



        %if k <= length(names)
        %    if k <= ci(1)
        %        % Associate specific name and color
        %        plot(x,u,'lineWidth',3,'DisplayName',names(k),'Color',colors(k,:))
        %    else
        %        % Associate specific name
        %        plot(x,u,'lineWidth',3,'DisplayName',names(k)) 
        %    end
        %elseif k <= ci(1)
        %    % Associate specific color
        %    temp = plot(x,u,'lineWidth',3,'Color',colors(k,:));
        %    set(get(get(temp, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        %else
        %    % Associate neither name nor color
        %    temp = plot(x,u,'lineWidth',3);
        %    set(get(get(temp, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        %end
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

    %axis([0 1 0 m.*1.1]);
    ax.XLim = [0 1];
    ax.YLim = [0 m*1.1];
    ylabel(ax,'Density');
    xlabel(ax,'Position');
    title(ax,tit,'Fontsize',18,'FontWeight', 'bold')
    %legend show

    return_data = ax;

end