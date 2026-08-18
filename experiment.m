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
                    'pr',1, ...     % pause rate
                    'ct',@(x) 1);   % cutoff function

% Add subdirectories to path
addpath(genpath(pwd));

disp('Ran Boilerplate')

%% Set up problem

p_width = 0.1;
p_eps = p_width/2;

parameters.s1    = p_eps;
parameters.s2    = 2*p_eps;
parameters.r1    = 1/3- p_eps;
parameters.r2    = 1/3 + 2*p_eps;
parameters.dt    = 0.001;
parameters.tfin  = 10;
parameters.del   = 0.1;
parameters.eps   = 0.015;
parameters.alph  = 0.9;
parameters.fr    = 20; 
parameters.pr    = 0.002;
parameters.ct    = ct_tanh_a(parameters.del);

disp('Set Parameters')

%% Generate Initial Conditions

% Setup domain
n = 15;
x = linspace(0,1,2^n);

c = 2*(1-4/3+6/5-4/7+1/9)/p_width;
build_pulse = @(x) @(z) (z>x-p_eps).*(z<x_p_eps).*(1-(z-x).^2).^4/c;

% Setup continuous IC: a k-peak wave with some thickness at the base
pulse_a = build_pulse(p_eps);
pulse_b = build_pulse(1/3+p_eps);
pulse_c = build_pulse(2/3+p_eps);

ic_s = pulse_a(x)+pulse_b(x)+pulse_c(x);
ic_d = ic_s/3;



disp('Generated IC')