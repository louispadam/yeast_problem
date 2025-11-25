%% 12/03/2025
%
% This script is for a research meeting with Keith Promislow.
%
% last updated 11/19/25 by Adam Petrucci

%% Boilerplate

addpath(genpath(pwd));
parameters = boiler_plate();

disp('Ran Boilerplate')

%% Set up problem

parameters.s1     = 0.1;
parameters.s2     = 0.3;
parameters.r1     = 0.6;
parameters.r2     = 0.7;
parameters.dt     = 0.001;
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
n = 9;
x = linspace(0,1,2^n);

% Setup initial: a k-peak wave with some thickness at the base
k = 4;
ic_s = @(z) (1-cos(2*pi*k*z))+0.5;
s_n = lp_integrate(x,ic_s(x),1);
ic = @(z) ic_s(z)/s_n;

ic_c = ic(x);
ic_p = sample(x,ic,1000);

disp("Generated Initial Data")

%% Visualize Initial Conditions

figure(1)
clf
ax = gca;
frame(x,parameters,ax,...
      "Data",{ic_c,ic_p},...
      "Meta",{struct('name',"Continuous IC", ...
                     'discrete',false, ...
                     'color',[0,0,255]/255) ...
              struct('name',"Particle IC", ...
                     'discrete',true, ....
                     'color',[100,100,255]/255, ...
                     'thickness',parameters.del)}, ...
      "Legend",true,...
      "Regions",true,...
      "Region_labels",true,...
      "Title","Initial Conditions");

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

figure(2)
clf
ax = gca;
frame(x,parameters,ax,...
      "Data",{fc_c,fc_p},...
      "Meta",{struct('name',"Continuous FC", ...
                     'discrete',false, ...
                     'color',[0,0,255]/255) ...
              struct('name',"Particle FC", ...
                     'discrete',true, ....
                     'color',[100,100,255]/255, ...
                     'thickness',parameters.del)}, ...
      "Legend",true,...
      "Regions",true,...
      "Region_labels",true,...
      "Title","Final Conditions");

disp('Visualized IC')

%% View Animation

figure(3)
clf

ax = gca;

animate(x,parameters,ax, ...
    "Data",{d_c,d_p},...
    "Meta",{struct('name',"Continuous FC", ...
                   'discrete',false, ...
                   'color',[0,0,255]/255) ...
            struct('name',"Particle FC", ...
                   'discrete',true, ....
                   'color',[100,100,255]/255, ...
                   'thickness',parameters.del)}, ...
    Time=t_c, ...
    Title="First Comparison", ...
    Legend=true, ...
    Regions=true, ...
    Region_labels=true);

disp("Finished Animating")