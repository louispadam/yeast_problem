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

num_particles = 1000;
ic_p = sort(sample(x,ic,num_particles));

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

%%

equi_data = metric_equilibrium(x,t_c,d_p,d_c);

disp('Computed equilibrium distance')

%%

opt_frame = figure(3);
clf(opt_frame);
ax = axes(opt_frame);

frame(x,parameters,ax,...
      "Data",{d_c(end-equi_data(2)+1,:),fc_p},...
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

sprintf('Minimum Wasserstein Distance of %0.4f',equi_data(1));

%%

i = 1;

parameters.s1 = s1_vals(i);

final_frame = figure(2);
clf(final_frame);
ax = axes(final_frame);

frame(x,parameters,ax,...
      "Data",{final_c(i,:),final_c(i,:)},...
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

disp('visualized')

%% Run Big Experiment

parameters.tfin   = 50;
parameters.update = false;

total_count = 0;

s1_vals = 0:0.05:0.25;
exp_colec = cell(1,length(s1_vals));
for s1 = 1:length(s1_vals)
    sl_vals = 0.1:0.05:0.25;
    exp_colec{s1} = cell(1,length(sl_vals));
    for sl = 1:length(sl_vals)
        rl_vals = 0.1:0.05:0.25;
        exp_colec{s1}{sl} = cell(1,length(rl_vals));
        for rl = 1:length(rl_vals)
            gp_vals = 0.1:0.05:(0.9-sl_vals(sl)-rl_vals(rl)-s1_vals(s1));
            exp_colec{s1}{sl}{rl} = cell(1,length(gp_vals));
            for gp = 1:length(gp_vals)
                exp_colec{s1}{sl}{rl}{gp} = cell([1,2]);
                exp_colec{s1}{sl}{rl}{gp}{1} = ...
                    [s1_vals(s1),...
                     s1_vals(s1) + sl_vals(sl),...
                     s1_vals(s1) + sl_vals(sl) + gp_vals(gp),...
                     s1_vals(s1) + sl_vals(sl) + + gp_vals(gp) + rl_vals(rl)];
                total_count = total_count + 1;

                parameters.s1 = exp_colec{s1}{sl}{rl}{gp}{1}(1);
                parameters.s2 = exp_colec{s1}{sl}{rl}{gp}{1}(2);
                parameters.r1 = exp_colec{s1}{sl}{rl}{gp}{1}(3);
                parameters.r2 = exp_colec{s1}{sl}{rl}{gp}{1}(4);

                [t_cont, d_cont] = cont_sudospec_etd_feuler(ic_c,parameters);
                [t_part, d_part] = particle_feuler(ic_p,parameters);

                exp_colec{s1}{sl}{rl}{gp}{2} = cell([1,2]);
                exp_colec{s1}{sl}{rl}{gp}{2}{1} = d_part(end,:);
                exp_colec{s1}{sl}{rl}{gp}{2}{2} = d_cont(end,:);

                update = sprintf("Finished running [ " + ...
                                 "s1 : %d/%d, " + ...
                                 "sl : %d/%d, " + ...
                                 "rl : %d/%d, " + ...
                                 "gp : %d/%d ]", s1,length(s1_vals),...
                                                 sl,length(sl_vals),...
                                                 rl,length(rl_vals), ...
                                                 gp,length(gp_vals));

                disp(update);

            end
        end
    end
end

disp('Completed experiment')

%%

save('big_exp_first_s1.mat','exp_colec','-v7.3')

disp('Done Saving')

%%

load("small_lol_exp.mat","exp_colec");

disp('Done Loading')

%%

exp_colec{1}{1}{1}{1}{1}

%%

total_count

%% Visualize Final Conditions

fc_c = exp_colec{1}{1}{1}{1}{2}{2};
fc_p = exp_colec{1}{1}{1}{1}{2}{1};

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

disp('Visualized FC')

%%

0:1:-3