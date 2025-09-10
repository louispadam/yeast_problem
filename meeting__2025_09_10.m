%% 09/10/2025
%
% This script is for a research meeting with Keith Promislow. It is
% intended to introduce the 'yeast problem' coding package and present
% early results.
%
% last updated 09/09/25 by Adam Petrucci

%% Boilerplate

% Instantiate Parameter Structure
parameters = struct('s1',0, ...     % beginning of signalling region
                    's2',0.3, ...   % end of signalling region
                    'r1',0.4, ...   % beginning of responsive region
                    'r2',0.5, ...   % end of responsive region
                    'dt',0.001, ... % simulation time step
                    'tfin',100, ... % ending time
                    'del',0.1, ...  % sharpness of hyperbolic cutoffs
                    'eps',0.025, ...% diffusive coefficient
                    'alph',0.9, ... % linear inhibition term
                    'fr',2, ...     % frame rate
                    'pr',1);        % pause rate

% Add subdirectories to path
addpath(genpath(pwd));

disp('Ran Boilerplate')

%% Set up problem

parameters.s1    = 0.1;
parameters.s2    = 0.3;
parameters.r1    = 0.6;
parameters.r2    = 0.7;
parameters.dt    = 0.001;
parameters.tfin  = 50;
parameters.del   = 0.1;
parameters.eps   = 0.025;
parameters.alph  = 0.9;
parameters.fr    = 20; 
parameters.pr    = 0.002;

disp('Set Parameters')

%% Generate Initial Conditions

% Setup domain
n = 11;
x = linspace(0,1,2^n);

% Setup continuous IC: a k-peak wave with some thickness at the base
k = 4;
ic_d = (1-cos(2*pi*k*x))+0.2;
ic_d = ic_d/lp_integrate(x,ic_d,1);

% Setup particle IC (by drawing from continuous IC)
trials = 3;
ic_p = zeros([trials,1000]);
for c = 1:trials
    ic_p(c,:) = construct_em(x,ic_d,1000);
end

disp('Generated IC')

%% Simulate Continuous System (MVP)

[ps_t, ps_d] = pseudospectral(ic_d,parameters);
ps_f = reshape(ps_d,[1,size(ps_d)]);

disp('Simulated MVP')

%% Simulate Particle System (NODE) with 250, 500, and 1000 particles
% Note to self: this section would be more convenient in a nested for loop

en_f_0250 = zeros([trials,length(ps_t),250]);
en_f_0500 = zeros([trials,length(ps_t),500]);
en_f_1000 = zeros([trials,length(ps_t),1000]);

for c = 1:trials

    p = ic_p(c,1:250);
    [en_t, en_p] = forward_euler_noise(p,parameters);
    en_f_0250(c,:,:) = reshape(en_p,[1,size(en_p)]);

    p = ic_p(c,1:500);
    [en_t, en_p] = forward_euler_noise(p,parameters);
    en_f_0500(c,:,:) = reshape(en_p,[1,size(en_p)]);

    p = ic_p(c,1:1000);
    [en_t, en_p] = forward_euler_noise(p,parameters);
    en_f_1000(c,:,:) = reshape(en_p,[1,size(en_p)]);

    disp("Ran trial " + c);

end

disp('Simulated NODE')

%% Graph Initial Conditions
% same note regarding for loop

th = zeros(1,trials);

for c = 1:trials
    m = metric_lp_1(en_f_0250(c,1,:),x,ps_f(1,1,:),1);
    th(c) = m(2);
end
frame(x,parameters,1,Continuous=ps_f,Discrete=en_f_0250,Thickness=th,Index=1);

for cc = 1:trials
    m = metric_lp_1(en_f_0500(c,1,:),x,ps_f(1,1,:),1);
    th(c) = m(2);
end
frame(x,parameters,2,Continuous=ps_f,Discrete=en_f_0500,Thickness=th,Index=1);

for cc = 1:trials
    m = metric_lp_1(en_f_1000(c,1,:),x,ps_f(1,1,:),1);
    th(c) = m(2);
end
frame(x,parameters,3,Continuous=ps_f,Discrete=en_f_1000,Thickness=th,Index=1);


disp('Graphed Initial Conditions')

%% Animate a Trial

animate(x,parameters,4,Continuous=ps_f,Discrete=en_f_1000(1,:,:),Thickness=[th(1)], ...
    Regions=true,Names=["ps","en1"]);

disp('Animated Trial 1')

%% Graph All Trials

frame(x,parameters,5,Continuous=ps_f,Discrete=en_f_1000,Thickness=th,Index=length(ps_f(1,:,1)))

disp('Graphed All Trials')

%%