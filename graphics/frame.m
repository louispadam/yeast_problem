function return_data = frame(x_data,parameters,fig_num,options)
%FRAME Animate a collection of simulations.
%
%last updated 08/31/25
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
    options.Time = 0                        % time point of simulation
    options.Index = 1                       % time step of simulation
    options.Continuous = []
    options.Discrete = []
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
    tt = options.Time;
    ind = options.Index;
    ps = options.Continuous;
    pp = options.Discrete;
    pt = options.Thickness;

    % Define Spatial Parameters
    a0 = params.s1;
    a1 = params.s2;
    b0 = params.r1;
    b1 = params.r2;
    del = params.del;

    %****************************
    % Compile Data
    %****************************
    
    ps_s = size(ps);
    pp_s = size(pp);

    y = zeros([pp_s(1)+ps_s(1),1,length(x)]);

    if ps_s(1) > 0
        y(1:ps_s(1),1,:) = ps(:,ind,:);
    end

    for k = 1:pp_s(1)
        y(ps_s(1)+k,1,:) = fatten_points_polynomial(x,reshape(pp(k,ind,:),[1,pp_s(3)]),pt(k));
    end

    %****************************
    % Construct Figures
    %****************************
    figure(fn)
    clf
    hold on

    si = size(y);
    ci = size(colors);
    for k = 1:si(1)
        u = reshape(y(k,1,:),[1,si(3)]);
        %u = squeeze(y(k,ind,:));
        if k <= length(names)
            if k <= ci(1)
                % Associate specific name and color
                plot(x,u,'lineWidth',3,'DisplayName',names(k),'Color',colors(k,:))
            else
                % Associate specific name
                plot(x,u,'lineWidth',3,'DisplayName',names(k)) 
            end
        elseif k <= ci(1)
            % Associate specific color
            temp = plot(x,u,'lineWidth',3,'Color',colors(k,:));
            set(get(get(temp, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        else
            % Associate neither name nor color
            temp = plot(x,u,'lineWidth',3);
            set(get(get(temp, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        end
    end
    
    % Determine upper bound of figure
    m = max(y(:,1,:),[],"all");

    ChiR = 1.1*m*0.25*(tanh(4*pi*(x-b0)/del)+1).*(tanh(4*pi*(b1-x)/del)+1);
    ChiS = 1.1*m*0.25*(tanh(4*pi*(x-a0)/del)+1).*(tanh(4*pi*(a1-x)/del)+1);

    if reg
        temp1 = plot(x,ChiR,'-','Color',[192,192,192]/255,'linewidth',3,'DisplayName','responsive');
        temp2 = plot(x,ChiS,'-','Color',[128,128,128]/255,'linewidth',3,'DisplayName','signaling');
        if ~regl
            set(get(get(temp1, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
            set(get(get(temp2, 'Annotation'), 'LegendInformation'),'IconDisplayStyle', 'off');
        end
    end

    axis([0 1 0 m.*1.1]);
    tlt=sprintf('Time =%2.5g',tt);
    ylabel('Density');
    xlabel('Position');
    title(tlt,'Fontsize',18)
    legend show

    return_data = 1;

end