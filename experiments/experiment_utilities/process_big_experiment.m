function return_data = process_big_experiment(data_f)
%PROCESS_BIG_EXPERIMENT presents results of big_experiment
%
%last updated 09/24/26 by Adam Petrucci

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Cluster Analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Analyze Clusters

    cluster_averages = mean(reshape(data_f{3}{1}(4,:), trials, []), 1);
    cluster_nn = reshape(cluster_averages, length(cluster_averages)/num_nums, []);

    mf_clustering = data_f{2}{1}(4,:);
    clustering_error = cluster_nn - mf_clustering.';

    num_part_seq = zeros([1,size(clustering_error,2)]);
    for j = 1:length(num_part_seq)
        num_part_seq(j) = length(data_f{3}{2}{j}(1,:));
    end

    disp('Finished Computing cluster averages')

    % Present Cluster Analysis

    cluster1_figure = figure(2);
    clf(cluster1_figure);
    ax = axes(cluster1_figure);

    hold on
    plot(ax,num_part_seq,clustering_error, ...
         'LineWidth',1);
    for k = -5:5
        plot(ax,num_part_seq,k*ones(size(num_part_seq)),...
             'Color',[100,100,100]/255, ...
             'LineStyle',':', ...
             'LineWidth',1)
    end
    hold off

    title(ax,"Average Cluster Error");
    xlabel(ax,"Number of Particles");
    ylabel(ax,"|NODE - MF|");
    ylim(ax,[-6,6]);

    disp('Presented Clustering-Error Figure')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Analyze Discrepancies Between NODE and MF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Compute error in final NOC

    [row_4, col] = find(abs(clustering_error) > 4);
    row_4 = unique(row_4);
    bad_sets_4 = data_f{1}(2:end,row_4);

    [row_3, col] = find(abs(clustering_error) > 3);
    row_3 = setdiff(unique(row_3),row_4);
    bad_sets_3 = data_f{1}(2:end,row_3);

    [row_2, col] = find(abs(clustering_error) > 2);
    row_2 = setdiff(unique(row_2),row_4);
    row_2 = setdiff(row_2,row_3);
    bad_sets_2 = data_f{1}(2:end,row_2);

    [row_1, col] = find(abs(clustering_error) > 1);
    row_1 = setdiff(unique(row_1),row_4);
    row_1 = setdiff(row_1,row_3);
    row_1 = setdiff(row_1,row_2);
    bad_sets_1 = data_f{1}(2:end,row_1);

    disp('Constructed index sets for NOC differences')

    % Collect full parameter set
    nonzero_params = data_f{1}(2:end,:);

    % Plot trends in error

    clustering2_figure = figure(3);
    clf(clustering2_figure);
    ax = axes(clustering2_figure);

    dark  = [065,105,225]/255;
    light = [137,207,240]/255;

    hold on
    scatter3(ax,bad_sets_4(1,:),bad_sets_4(2,:),bad_sets_4(3,:), ...
             'filled','MarkerFaceColor',dark*4/4+light*0/4, ...
             'DisplayName','Error >4');
    scatter3(ax,bad_sets_3(1,:),bad_sets_3(2,:),bad_sets_3(3,:), ...
             'filled','MarkerFaceColor',dark*3/4+light*1/4, ...
             'DisplayName','Error >3');
    scatter3(ax,bad_sets_2(1,:),bad_sets_2(2,:),bad_sets_2(3,:), ...
             'filled','MarkerFaceColor',dark*2/4+light*2/4, ...
             'DisplayName','Error >2');
    scatter3(ax,bad_sets_1(1,:),bad_sets_1(2,:),bad_sets_1(3,:), ...
             'filled','MarkerFaceColor',dark*1/4+light*3/4, ...
             'DisplayName','Error >1');
    scatter3(ax,nonzero_params(1,:),nonzero_params(2,:),nonzero_params(3,:),...
             'MarkerEdgeColor','black', ...
             'DisplayName','Parameter Space');
    hold off

    title(ax,'Trends in Error')
    legend(ax)

    xlabel(ax,'s_2')
    ylabel(ax,'r_1')
    zlabel(ax,'r_2')

    xlim(ax,[0, 1])
    ylim(ax,[0, 1])
    zlim(ax,[0, 1])

    disp('Presented Trends in Clustering Error')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Analyze trends in convergence
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Check convergence trends

    no_conv_mf = find(~data_f{2}{1}(1,:));
    no_conv_mf_params = data_f{1}(2:end,no_conv_mf);

    disp('Collected Convergence Data')

    % Present convergence data

    convergence1_figure = figure(4);
    clf(convergence1_figure);
    ax = axes(convergence1_figure);

    hold on
    scatter3(ax,no_conv_mf_params(1,:), ...
                no_conv_mf_params(2,:), ...
                no_conv_mf_params(3,:), ...
                'filled','MarkerFaceColor','black', ...
                'DisplayName','No Convergence');
    scatter3(ax,nonzero_params(1,:),nonzero_params(2,:),nonzero_params(3,:),...
             'MarkerEdgeColor','black', ...
             'DisplayName','Parameter Space');
    hold off

    title(ax,'Trends in Convergence')
    legend(ax)

    xlabel(ax,'s_2')
    ylabel(ax,'r_1')
    zlabel(ax,'r_2')

    xlim(ax,[0, 1])
    ylim(ax,[0, 1])
    zlim(ax,[0, 1])

    disp('Presented Trends in Convergence')

    return_data = 1;

end