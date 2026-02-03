%% 01/20/2026
%
% This script is for a research meeting with Keith Promislow.
%
% last updated 01/17/26 by Adam Petrucci

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
parameters.tfin   = 50;
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

disp('Simulated MVP')

%% Simulate Particle System (NODE)

[t_p, d_p] = particle_feuler(ic_p,parameters);

disp('Simulate NODE')

%% Visualize Final Conditions

fc_c = d_c(end,:);
fc_p = d_p(end,:);

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

animation_fig = figure(3);
clf(animation_fig);
ax = axes(animation_fig);

animate(x,parameters,ax, ...
    Data={d_c,d_p}, ...
    Meta={struct('name','Continuous', ...
                'discrete',false, ...
                'color',[253,179,102]/255) ...
          struct('name',"Particle", ...
                'discrete',true, ....
                'color',[110,166,205]/255, ...
                'thickness',parameters.del)}, ...
    Time=t_p, ...
    Title="Time Series", ...
    Legend=true, ...
    Regions=true, ...
    Region_labels=true);

disp("Finished Animating")

%%

d_p_d = histcounts(fc_p, x);
d_p_d = cat(2,d_p_d,[0]);
lp = lp_integrate(x,d_p_d,1);
d_p_d = d_p_d/lp;

steps_in_rev = sum(t_c>t_c(end)-1);

trans_dists = zeros([1,steps_in_rev]);

for i = 1:steps_in_rev
     trans_dists(i) = lp_integrate(x,abs(d_c(end-i+1,:)-d_p_d),1);
end

[min_val, min_ind] = min(trans_dists);

%%

figure(11)
plot(1:steps_in_rev,trans_dists)

%%

opt_frame = figure(3);
clf(opt_frame);
ax = axes(opt_frame);

frame(x,parameters,ax,...
      "Data",{d_c(end-min_ind+1,:),fc_p},...
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
      "Title","Optimal Overlap");

sprintf('Minimum Wasserstein Distance of %0.4f',min_val);