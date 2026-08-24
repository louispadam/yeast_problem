function return_data = present_cluster_detection(time,data,options)
arguments (Input)
    time
    data
end
arguments (Input)
    options.Figure = 1
    options.Title = 'Clustering Modes'
end
arguments (Output)
    return_data
end

    to_track = size(data,2);

    clustering_frame = figure(options.Figure);
    clf(clustering_frame);
    ax = axes(clustering_frame);

    c = turbo(16);
    blueColors = c(1:to_track,:);
    ind_freq = round(linspace(1,length(time),10));

    hold on
    for i = 1:to_track
        plot(ax,time,data(:,i), ...
                LineWidth=2, ...
                Color=blueColors(i,:),...
                DisplayName=sprintf('Mode = %d',i));
        text(time(ind_freq), data(ind_freq,i),...
            string(i*ones([1,length(ind_freq)])), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle',...
            'Color',blueColors(i,:),...
            'FontSize',8,...
            'BackgroundColor','white', ...
            'Margin',0.5)
    end
    hold off
    legend(Location="west",NumColumns=1)

    title(ax,options.Title)
    xlabel(ax,'Time')
    ylabel(ax,'Magnitude of Mode')

    return_data = 1;

end