%% Boilerplate

% Add project to path
addpath(genpath(pwd));

% Build parameter structures
parameters = instantiate_parameters();

disp('Ran Boilerplate')

%% Set up problem

parameters.s1     = 0.10;
parameters.s2     = 0.30;
parameters.r1     = 0.65;
parameters.r2     = 0.80;
parameters.eps    = 0;
parameters.alph   = 0.9;
parameters.ct     = ct_sharp(0.01);

disp('Set Parameters')

%% Generate Initial Data

% Setup domain
x = linspace(0,1-2^(-16),2^16);

k = 4;
ic_s = @(z) (1-cos(2*pi*k*z))+0.5;
ic_n = ic_s(x);
s_n = 0.5 * sum([(ic_n(2:end) + ic_n(1:end-1)) .* (x(2:end) - x(1:end-1)), ...
                 (ic_n(1) + ic_n(end)) * mod(x(1) - x(end),1)]);
ic = @(z) ic_s(z)/s_n;

generate_sample = construct_sampler(ic);
ic_c = ic(x);

lag_data = euler_to_lagrange(ic_mesh,ic_c);


disp('Prepared Initial Data')

%% For Research Statement

% Run mean-field simulation
[time_mf, mesh_mf, mass_mf, ~, conv_mf, end_mf] = ...
                         cont_cons_lag(lag_data,parameters,...
                                       "EndTime",Inf,...
                                       "Collect",true,...
                                       "Track",true);

%%

fc_mass = circshift(ic_mass,-mass_mf(end));
fc_eul = lagrange_to_euler(squeeze(mesh_mf(end,:)),x,fc_mass);

disp('done')

%%

ic_p = generate_sample(1000);

[time_node, data_node, ~, ~, end_node] = ...
                            particle_proof(ic_p,parameters, ...
                                           "Update",true,...
                                           "Collect",true);

%%

node_color = [31,119,180]/255;
mefi_color = [255,127,14]/255;

ic_figure = figure(1);
clf(ic_figure);
ic_ax = gca(ic_figure);

ic_figure.Units = 'inches';
ic_figure.Position = [1 1 5.5 3.5];

frame(x,parameters,ic_ax,Data={ic_p,...
                               ic_c},...
                              Meta={struct(name="NODE",...
                                           discrete=true,...
                                           color=node_color,...
                                           thickness=0.01),...
                                    struct(name="MF",...
                                           discrete=false,...
                                           color=mefi_color)},...
                              Legend=true,...
                              Regions=true,...
                              Region_labels=true,...
                              Title="Initial State");

%exportgraphics(ic_figure,'ic_figure.pdf','ContentType','vector')

fc_figure = figure(2);
clf(fc_figure);
fc_ax = gca(fc_figure);

fc_figure.Units = 'inches';
fc_figure.Position = [1 1 5.5 3.5];

frame(x,parameters,fc_ax,Data={squeeze(data_node(end,:)),...
                               fc_eul},...
                              Meta={struct(name=sprintf('NODE (t=%2.0f)',time_node(end)),...
                                           discrete=true,...
                                           color=node_color,...
                                           thickness=0.000013),...
                                    struct(name=sprintf('MF (t=%2.0f)',time_mf(end)),...
                                           discrete=false,...
                                           color=mefi_color)},...
                              Legend=true,...
                              Regions=true,...
                              Region_labels=true,...
                              Title="Final State");

%exportgraphics(fc_figure,'fc_figure.pdf','ContentType','vector')

%% Little experiment

modes = 7;

% Run mean-field simulation
[time_mf, data_mf, ~, ~, end_mf] = ...
                         cont_cons_lag(ic_c,parameters,...
                                       "Collect",true);

% Detect clusters of mean-field simulation
clusters_mf = track_clusters(time_mf,data_mf,"Modes",modes,"Discrete",false);

% Derive cluster data
[trans_mf,stable_mf] = analyze_clusters(time_mf,clusters_mf);

%%

trials = 15;
N = 1000*(1:20);
exp_colec = cell([size(N),1]);

for n = 1:length(N)

    num_par = N(n);

    exp_colec{n} = struct('num_par',num_par,...
                          'node_endstate',zeros([trials,1]),...
                          'node_endtime',zeros([trials,1]),...
                          'node_enddata',zeros([trials,N(n)]),...
                          'node_transtime',zeros([trials,1]),...
                          'node_clusters',zeros([trials,1]));

    for t = 1:trials

        % Sample initial particle data
        ic_p = sort(generate_sample(num_par));

        % Run NODE simulation
        [time_node, data_node, ~, ~, end_node] = ...
                            particle_proof(ic_p,parameters, ...
                                           "Update",true,...
                                           "Collect",true);

        % Detect clusters of NODE simulation
        clusters_node = track_clusters(time_node,data_node,...
                                               "Modes",modes,...
                                               "Discrete",true);

        % Derive cluster data
        [trans_node,stable_node] = analyze_clusters(time_node,clusters_node);

        exp_colec{n}.node_endstate(t) = end_node;
        exp_colec{n}.node_endtime(t) = time_node(end);
        exp_colec{n}.node_enddata(t,:) = squeeze(data_node(end,:));
        exp_colec{n}.node_transtime(t) = trans_node;
        exp_colec{n}.node_clusters(t) = stable_node;

    end % end of NODE trials

end % for N

save("26_10_08_experiment.mat","exp_colec");

disp('done NODE')

%%

all_transtime = zeros([length(N),trials]);
all_endtime = zeros([length(N),trials]);
for n = 1:length(N)
    all_transtime(n,:) = exp_colec{n}.node_transtime;
    all_endtime(n,:) = exp_colec{n}.node_endtime;
end

transtime_figure = figure(5);
clf(transtime_figure);
transtime_ax = gca(transtime_figure);

hold on
for n = 1:length(N)
    scatter(transtime_ax,N(n)*ones([1,trials]),all_transtime(n,:));
end
plot(transtime_ax,N,trans_mf*ones([length(N),1]));
title(transtime_ax,'Transtime');

transtime_figure.Units = 'inches';
transtime_figure.Position = [1 1 6.0 2.5];

%%

endtime_figure = figure(6);
clf(endtime_figure);
endtime_ax = gca(endtime_figure);

process_endtime = sort(all_endtime,2);
process_endtime = process_endtime(:,1:10);
avg_endtime = mean(process_endtime,2);

hold on
plot_all = scatter(endtime_ax,N,process_endtime,'k');
plot_avg = scatter(endtime_ax,N,avg_endtime,'k','filled');
plot_mf = plot(endtime_ax,N,time_mf(end)*ones([length(N),1]),'--r');

title(endtime_ax,'Convergence Time Comparison of NODE and MF');
legend([plot_mf,plot_avg], {sprintf('MF Convergence Time (t=%3.0f)',time_mf(end)),'NODE Average Convergence Time'})
xlabel('Number of Particles')
ylabel('Convergence Time')
xlim([0 N(end)*1.1])
ylim([time_mf(end)*0.3 time_mf(end)*1.2])
endtime_ax.XAxis.Exponent = 0;
xticks(0:5000:20000)
xtickformat('%.0f')
yticks(50:25:200)
ytickformat('%.0f')

saveas(endtime_figure, 'endtime.png')
endtime_figure.Units = 'inches';
endtime_figure.Position = [1 1 6.0 2.5];

disp('done endtime plot')


%% Running some code

node_color = [1,0,0];
mefi_color = [0,1,0];

num_par = 1000;

ic_p = sort(generate_sample(num_par));

[time_node, data_node] = particle_proof(ic_p,...
                                        parameters,...
                                        "Update",true,...
                                        "EndTime",Inf);
final_node = sort(squeeze(data_node(end,:)));

[time_mf, data_mf, ~, conv_mf] = cont_cons_lag(ic_c,...
                                   parameters,...
                                   "Update",true,...
                                   "Collect",true, ...
                                   "EndTime",Inf, ...
                                   "Track",true);
final_lag_mefi = squeeze(data_mf(end,:,:));
final_eul_mefi = lagrange_to_euler(squeeze(final_lag_mefi(1,:)),x,...
                                   squeeze(final_lag_mefi(2,:)));

disp('finished simulations')

%%

conv_figure = figure(5);
clf(conv_figure);
conv_ax = gca(conv_figure);

plot(conv_ax,time_mf(2:end),conv_mf)

%% Making some figures

ic_lag_figure = figure(1);
clf(ic_lag_figure);
ic_lag_ax = gca(ic_lag_figure);

frame_cumsum(parameters,ic_lag_ax,Data={[ic_p;ones(size(ic_p))/num_par],...
                                        [x;ic_c.*([x(2:end),x(1)+1]-x)]},...
                              Meta={struct(name="node",...
                                           color=node_color),...
                                    struct(name="mefi",...
                                           color=mefi_color)},...
                              Legend=true,...
                              Regions=true,...
                              Region_labels=true,...
                              Title="IC");

fc_lag_figure = figure(2);
clf(fc_lag_figure);
fc_lag_ax = gca(fc_lag_figure);

frame_cumsum(parameters,fc_lag_ax,Data={[final_node;ones(size(final_node))/num_par],...
                                     final_lag_mefi},...
                              Meta={struct(name="node",...
                                           color=node_color),...
                                    struct(name="mefi",...
                                           color=mefi_color)},...
                              Legend=true,...
                              Regions=true,...
                              Region_labels=true,...
                              Title="FC lag");

fc_eul_figure = figure(3);
clf(fc_eul_figure);
fc_eul_ax = gca(fc_eul_figure);

frame(x,parameters,fc_eul_ax,Data={final_node,...
                                 final_eul_mefi},...
                              Meta={struct(name="node",...
                                           discrete=true,...
                                           color=node_color,...
                                           thickness=0.00001),...
                                    struct(name="mefi",...
                                           discrete=false,...
                                           color=mefi_color)},...
                              Legend=true,...
                              Regions=true,...
                              Region_labels=true,...
                              Title="FC eul");

disp('finished figures')

%%

make_video(x,parameters,figure(4),...
           "Data",{data_mf},...
           "Meta",{struct(name = "Check",...
                          discrete = false,...
                          color = [0,0,0],...
                          thickness = 0.01)},...
           "Time",time_mf,...
           "Regions",true,...
           "Region_labels",true,...
           "Legend",true);

%% Load in data

clear exp_data exp_error_node exp_error_mf

exp_data = load("data_public\26_08_27\big_experiment_26_08_27.mat");
exp_data = exp_data.exp_colec;
exp_error_node = load("data_public\26_08_27\error_node_big_experiment_26_08_27.mat");
exp_error_node = exp_error_node.error_data_node;
exp_error_mf = load("data_public\26_08_27\error_mf_big_experiment_26_08_27.mat");
exp_error_mf = exp_error_mf.error_data_mf;

disp('Done Loading Data')

%%

for e = 1:length(exp_error_mf)
    err_params = exp_error_mf(e).parameters;
    fprintf("Error with parameter set\n" + ...
            "[%.2f, %.2f, %.2f, %.2f]\n", ...
            err_params.s1, err_params.s2, err_params.r1, err_params.r2);
end

disp('Processed Errors')

%%

for ind1 = 1:length(exp_data)
    for ind2 = 1:length(exp_data{ind1})
        for ind3 = 1:length(exp_data{ind2})
            try
                to_print = exp_data{ind1}{ind2}{ind3}{1};
                fprintf("s1 = %.2f   s2 = %.2f   r1 = %.2f   r2 = %.2f\n", ...
                         to_print.s1, to_print.s2, to_print.r1, to_print.r2);
            catch
                fprintf("Vector error\n");
            end
        end
    end
end

%% Reorganize data structure
% Figure out another time

% Might be explicit formula for this, I should look into it later
total_sets = 0;
for ind1 = 1:length(exp_data)
    for ind2 = 1:length(exp_data{ind1})
        for ind3 = 1:length(exp_data{ind1}{ind2})
            total_sets = total_sets + 1;
        end
    end
end

num_nums = length(exp_data{1}{1}{1}{3});
trials = length(exp_data{1}{1}{1}{3}{1});

data_f = cell([1,3]);
data_f{2} = cell([1,2]);
data_f{3} = cell([1,2]);
data_f{3}{1} = NaN([4,trials * total_sets * num_nums]);
data_f{3}{2} = cell([1,num_nums]);

for k = 1:num_nums
    data_f{3}{2}{k} = NaN([trials,length(exp_data{1}{1}{1}{3}{k}{1}{3})]);
end

param_array = NaN([4,total_sets]);
mf_stats = NaN([4,total_sets]);
mf_states = NaN([length(x),total_sets]);

temp_check = 1;

cc = 1;
for ind1 = 1:length(exp_data)
    for ind2 = 1:length(exp_data{ind1})
        for ind3 = 1:length(exp_data{ind1}{ind2})

            if isempty(exp_data{ind1}{ind2}{ind3}{2}{4})
                mf_stat = [exp_data{ind1}{ind2}{ind3}{2}{1},...
                           exp_data{ind1}{ind2}{ind3}{2}{2},...
                           -1,...
                           exp_data{ind1}{ind2}{ind3}{2}{5}];
                temp_check = temp_check + 1;
                if temp_check == 3
                    ind1
                    ind2
                    ind3
                end
            else
                mf_stat = [exp_data{ind1}{ind2}{ind3}{2}{1},...
                           exp_data{ind1}{ind2}{ind3}{2}{2},...
                           exp_data{ind1}{ind2}{ind3}{2}{4},...
                           exp_data{ind1}{ind2}{ind3}{2}{5}];
            end

            mf_state = exp_data{ind1}{ind2}{ind3}{2}{3};

            param_set = exp_data{ind1}{ind2}{ind3}{1};
            param_array(:,cc) = [param_set.s1,...
                                 param_set.s2,...
                                 param_set.r1,...
                                 param_set.r2];
            mf_stats(:,cc) = mf_stat;
            mf_states(:,cc) = mf_state;

            for n = 1:num_nums
                for t = 1:trials

                    data_f{3}{2}{n}(t,:) = exp_data{1}{1}{1}{3}{n}{t}{3};

                    if isempty(exp_data{ind1}{ind2}{ind3}{3}{n}{t}{4})
                        node_stat = [exp_data{ind1}{ind2}{ind3}{3}{n}{t}{1},...
                                 exp_data{ind1}{ind2}{ind3}{3}{n}{t}{2},...
                                 -1,...
                                 exp_data{ind1}{ind2}{ind3}{3}{n}{t}{5}];
                    else
                        node_stat = [exp_data{ind1}{ind2}{ind3}{3}{n}{t}{1},...
                                     exp_data{ind1}{ind2}{ind3}{3}{n}{t}{2},...
                                     exp_data{ind1}{ind2}{ind3}{3}{n}{t}{4},...
                                     exp_data{ind1}{ind2}{ind3}{3}{n}{t}{5}];
                    end

                    data_f{3}{1}(:,(n-1)*total_sets*trials + ...
                                   (cc-1)*trials + t) = node_stat;

                end
            end

            cc = cc + 1;
        end
    end
end

data_f{1} = param_array;
data_f{2}{1} = mf_stats;
data_f{2}{2} = mf_states;

disp('Done reorganizing data')

%%

cluster_averages = mean(reshape(data_f{3}{1}(4,:), trials, []), 1);
cluster_nn = reshape(cluster_averages, length(cluster_averages)/num_nums, []);

mf_clustering = data_f{2}{1}(4,:);
clustering_error = cluster_nn - mf_clustering.';

num_part_seq = zeros([1,size(clustering_error,2)]);
for j = 1:length(num_part_seq)
    num_part_seq(j) = length(data_f{3}{2}{j}(1,:));
end

disp('Finished Computing cluster averages')

%% Cluster Analysis Presentation

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

%% Analyze trends in error

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

disp('Constructed index sets')

%% Collect full parameter set

nonzero_params = data_f{1}(2:end,:);

%% Plot trends in error

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

%% Check convergence trends

no_conv_mf = find(~data_f{2}{1}(1,:));
no_conv_mf_params = data_f{1}(2:end,no_conv_mf);

disp('Collected Convergence Data')

%% Present convergence data

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