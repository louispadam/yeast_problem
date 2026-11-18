function return_data = analyze_node_data(mf_data,exp)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    mf_transtime = mf_data(1);
    mf_endtime = mf_data(2);
    mf_noc = mf_data(3);

    trials = length(exp{1}.node_data{1}.node_endstate);
    N = zeros([1,length(exp{1}.node_ic)]);
    for n = 1:length(exp{1}.node_ic)
        N(n) = exp{1}.node_data{n}.num_par;
    end

    % NODE data

    % Initalize storage objects
    compile_transtime = zeros(length(N),trials);
    compile_endtime = zeros(length(N),trials);
    compile_noc = zeros(length(N),trials);

    % Retrieve data from experiment
    for n = 1:length(N)
        for t = 1:trials
            compile_transtime(n,:) = exp{1}.node_data{n}.node_transtime;
            compile_endtime(n,:) = exp{1}.node_data{n}.node_endtime;
            compile_noc(n,:) = exp{1}.node_data{n}.node_clusters;
        end
    end

    % Perform statistical analysis
    average_transtime = mean(compile_transtime);
    average_endtime = mean(compile_endtime);
    agreement = sum(compile_noc == mf_noc) / trials;

    % Create figure for transition time

    transtime_figure = figure(1);
    clf(transtime_figure);
    transtime_ax = axes(transtime_figure);

    hold on
    scatter(transtime_ax,N,average_transtime,"black","filled",Marker="o")
    scatter(transtime_ax,N,compile_transtime.',"blue",Marker="o")
    plot(transtime_ax,[0,N(end)*1.1],mf_transtime*ones([1,2]),'--r');
    hold off
    xlim(transtime_ax,[0,N(end)*1.1])
    ylim(transtime_ax,[min([mf_transtime,average_transtime])*0.9,...
                       max([mf_transtime,average_transtime])*1.1])

    % Create figure for end time

    endtime_figure = figure(2);
    clf(endtime_figure);
    endtime_ax = axes(endtime_figure);

    hold on
    scatter(endtime_ax,N,average_endtime,"black","filled",Marker="o")
    scatter(endtime_ax,N,compile_endtime.',"blue",Marker="o")
    plot(endtime_ax,[0,N(end)*1.1],mf_endtime*ones([1,2]),'--r');
    hold off
    xlim(endtime_ax,[0,N(end)*1.1])
    ylim(endtime_ax,[min([mf_endtime,average_endtime])*0.9,...
                     max([mf_endtime,average_endtime])*1.1])

    % Create figure for NOC

    noc_figure = figure(3);
    clf(noc_figure);
    noc_ax = axes(noc_figure);

    hold on
    plot(noc_ax,N,agreement,'-ob');
    hold off

    return_data = 1;

end