%% 12/03/2025
%
% This script is for a research meeting with Keith Promislow. It is
% intended to compare results of simulating continuum and particle models
% now that the issue of numerical diffusion has been effectively resolved.
%
% last updated 11/30/25 by Adam Petrucci

%% Boilerplate

addpath(genpath(pwd));
parameters = boiler_plate();

disp('Ran Boilerplate')

%% Set up problem

parameters.s1     = 0.1;
parameters.s2     = 0.3;
parameters.r1     = 0.6;
parameters.r2     = 0.7;
parameters.dt     = 0.0004;
parameters.tfin   = 100;
parameters.del    = 0.01;
parameters.eps    = 0.0;      % NO DIFFUSION
parameters.alph   = 0.9;
parameters.fr     = 20; 
parameters.pr     = 0.002;
parameters.ct     = ct_tanh_k(parameters.del);
parameters.m_sz   = 2^15;
parameters.update = true;

disp('Set Parameters')

%% Generate Initial Data

% Setup domain
n = 11;
x = linspace(0,1,2^n);

% Setup initial: a k-peak wave with some thickness at the base
k = 4;
ic_s = @(z) (1-cos(2*pi*k*z))+0.5;
s_n = lp_integrate(x,ic_s(x),1);
ic = @(z) ic_s(z)/s_n;

ic_c = ic(x);
ic_p = sort(sample(x,ic,1000));

disp("Generated Initial Data")

%% Visualize Initial Conditions

initial_frame = figure(1);
clf(initial_frame);
ax = axes(initial_frame);

frame(x,parameters,ax,...
      "Data",{ic_c,ic_p},...
      "Meta",{struct('name',"Continuous IC", ...
                     'discrete',false, ...
                     'color',[253,179,102]/255) ...
              struct('name',"Particle IC", ...
                     'discrete',true, ....
                     'color',[110,166,205]/255, ...
                     'thickness',parameters.del)}, ...
      "Legend",true,...
      "Regions",true,...
      "Region_labels",true,...
      "Title","Initial Conditions");

%saveas(initial_frame,'intial_data_4peak.png')

disp('Visualized IC')

%% Simulate Continuous System (MVP)

[t_c, d_c] = cont_sudospec_etd_feuler(ic_c,parameters);

fc_c = d_c(end,:);

disp('Simulated MVP')

%% Simulate Particle System (NODE)

[t_p, d_p] = particle_feuler(ic_p,parameters);

fc_p = d_p(end,:);

disp('Simulate NODE')

%% Visualize Final Conditions

final_frame = figure(2);
clf(final_frame);
ax = axes(final_frame);

frame(x,parameters,ax,...
      "Data",{fc_c,fc_p},...
      "Meta",{struct('name',"Continuous FC", ...
                     'discrete',false, ...
                     'color',[253,179,102]/255) ...
              struct('name',"Particle FC", ...
                     'discrete',true, ....
                     'color',[110,166,205]/255, ...
                     'thickness',parameters.del)}, ...
      "Legend",true,...
      "Regions",true,...
      "Region_labels",true,...
      "Title","Final Conditions");

%saveas(final_frame,'final_data_4peak.png')

disp('Visualized IC')

%% View Animation

% animation_fig = figure(3);
% clf(animation_fig);
% ax = axes(animation_fig);
% 
% animate(x,parameters,ax, ...
%     "Data",{d_c,d_p},...
%     "Meta",{struct('name',"Continuous FC", ...
%                    'discrete',false, ....
%                    'color',[0,0,255]/255), ...
%             struct('name',"Particle FC", ...
%                    'discrete',true, ....
%                    'color',[100,100,255]/255, ...
%                    'thickness',parameters.del)}, ...
%     Time=t_c, ...
%     Title="First Comparison", ...
%     Legend=true, ...
%     Regions=true, ...
%     Region_labels=true);
% 
% disp("Finished Animating")

%% Measure Convergence

freq = 50;
keep = floor(1/(freq*parameters.dt));
if keep == 0
    keep = 1;
end

tot = floor(length(t_c) / keep);
was = zeros([1,tot]);

for k = 1:tot 
    
    j = keep*k;

    d_p_d = histcounts(d_p(j,:), x);
    lp = lp_integrate(x,d_p_d,1);
    d_p_d = d_p_d/lp;

    w = metric_wasserstein_new2(x,d_c(k,:),d_p_d);

    was(k) = w(1);

    fprintf("Calculated Iterate %d / %d\n",j,length(t_c));

end

%% Visualize Wasserstein Convergence

was_dist = figure(4);
clf('reset',was_dist)
plot((keep*(1:tot)+1)/length(t_c)*parameters.tfin,was)
title("Wasserstein Distance Between Continuum and Particle Models");
ylabel('Distance');
xlabel('Time');

%saveas(was_dist,'was_dist_4peak.png')

%% Consider Long-term (Equilibrium) Behavior

equi = d_p(floor(length(d_p)/2):end,:);
equi_d_1 = dist_per(equi(1,1),equi(1,2:end),1);

%% Visualize Long-term relative particle distances in space in NODE

per_dist1 = figure(5);
clf('reset',per_dist1)
plot(2:length(equi(1,:)),equi_d_1)
title("Periodic Distance in Particle System");
ylabel('Distance');
xlabel('Time');

%saveas(per_dist1,'per_dist1.png')

%%

d_c_cum = zeros([1,length(x)]);
for k = 2:length(fc_c)
    d_c_cum(k) = lp_integrate(x(1:k),fc_c(1:k),1);
end
d_c_cum(1) = 0;
d_c_cum(end) = 1;

d_c_cum_fig= figure(6);
clf('reset',d_c_cum_fig)
plot(x,d_c_cum)
title("Cumulative Distribution of Continuous System");
ylabel('Distance');
xlabel('Time');

%saveas(d_c_cum_fig,'d_c_cum.png')

%% Compute Long-term relative particle distances in time

p0 = equi(:,1);
p1 = equi(:,200);
p2 = equi(:,600);

compare01 = dist_per(p0,p1,1);
compare02 = dist_per(p0,p2,1);

%% Visualize long-term relative particle distances in time

per_analy_p = figure(7);
clf('reset',per_analy_p)
hold on

plot(1:length(equi(end-1000:end,1)),compare01(end-1000:end))
plot(1:length(equi(end-1000:end,1)),compare02(end-1000:end))
legend("B","C")
title("Periodic Analysis of Particles")
ylabel("Periodic Distance")
xlabel("Time")

%saveas(per_analy_p,'per_analy_p.png')

%%

per_analy_c = figure(8);
clf(per_analy_c,'reset')
hold on

hang_on = 0;

all_peaks = zeros([1000,3]);
for k = 1:1000
    current_step = d_c(length(t_c)-1000+k,:);
    use_step = cat(2,cat(2,current_step(end-100:end-1),current_step),current_step(2:100));
    corr_x = cat(2,cat(2,x(end-100:end-1)-1,x),x(2:100)+1);
    [pks, loc] = findpeaks(use_step,corr_x,'MinPeakHeight',5,...
                                          'MinPeakDistance',0.1);
    corr_loc = loc(0 <= loc & loc < 1);
    all_peaks(k,:) = corr_loc;
end

all_peaks_corr = all_peaks;
for k = 2:1000
    if all_peaks(k,1) < 0.1 && all_peaks(k-1,1) > 0.1
        all_peaks_corr(k:end,:) = circshift(all_peaks_corr(k:end,:),-1,2);
    end
end

all_dists = zeros([1000,2]);
for k = 1:1000
    loc = all_peaks_corr(k,:);
    for j = 2:length(corr_loc)
        all_dists(k,j-1) = dist_per(loc(1),loc(j),1);
        all_dists(j-1) = dist_per(corr_loc(1),corr_loc(j),1);
    end
end


plot(all_dists)
legend("B","C")
title("Periodic Analysis of Continuum")
ylabel("Periodic Distance")
xlabel("Time")

%saveas(per_analy_c,'per_analy_c.png')

%% Figure for Website

website_frame = figure(20);
clf(website_frame);
ax = axes(website_frame);

% [253,179,102]/255
frame(x,parameters,ax,...
      "Data",{ic_p,d_p(5000,:)},...
      "Meta",{struct('name',"Initial State", ...
                     'discrete',true, ...
                     'color',[253,149,72]/255, ...
                     'thickness',parameters.del)...
              struct('name',"Ending State", ...
                     'discrete',true, ....
                     'color',[110,166,205]/255, ...
                     'thickness',parameters.del)}, ...
      "Legend",true,...
      "Regions",false,...
      "Region_labels",false,...
      "Title","Particle Clustering in Yeast");

p = ax.Position;
a = annotation('textbox', ...
    [p(1)+0.15*p(3),p(2)+0.8*p(4), 0.1, 0.1], ...
    'String', sprintf('Transition:\n4-peak  >>  3-peak'), ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontWeight', 'bold', ...
    'FontSize', 11);

%saveas(website_frame,'graphic_yeast.png')

disp('for site')