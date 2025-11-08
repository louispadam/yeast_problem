%% Boilerplate

% Instantiate Parameter Structure
parameters = struct('s1',0, ...        % beginning of signalling region
                    's2',0.3, ...      % end of signalling region
                    'r1',0.4, ...      % beginning of responsive region
                    'r2',0.5, ...      % end of responsive region
                    'dt',0.001, ...    % simulation time step
                    'tfin',100, ...    % ending time
                    'del',0.1, ...     % sharpness of hyperbolic cutoffs
                    'eps',0.025, ...   % diffusive coefficient
                    'alph',0.9, ...    % linear inhibition term
                    'fr',2, ...        % frame rate
                    'pr',1, ...        % pause rate
                    'ct',@(x) 1, ...   % cutoff function
                    'm_sz',1);         % max vector size

% Add subdirectories to path
addpath(genpath(pwd));

disp('Ran Boilerplate')

%% Set up problem

parameters.s1    = 0.45;      % very small and separated regions should 
parameters.s2    = 0.5;       % have little impact on mechanics, in
parameters.r1    = 0.1;       % particular diffusion
parameters.r2    = 0.15;
parameters.dt    = 0.001;
parameters.tfin  = 10;
parameters.del   = 0.01;
parameters.eps   = 0.00;      % NO DIFFUSION
parameters.alph  = 0.00;      % NO INTERACTION
parameters.fr    = 1; 
parameters.pr    = 0.002;
parameters.ct    = ct_tanh_a(parameters.del);
parameters.m_sz  = 2^15;

disp('Set Parameters')

%% Generate Initial Conditions

% Setup domain
n = 11;
x = linspace(0,1,2^n);

% Set up initial distrubution: a narrow gaussian
start_width = 0.005;
ic_s = @(z) exp(-(z-0.5).^2/(start_width^2));
lp = lp_integrate(x,ic_s(x),1);
ic = @(z) ic_s(z)/lp;

ic_p = ic(x);

disp("Generated Initial Conditions")

%% Visualize Initial Conditions

figure(1)
clf
ax = gca;
frame(x,parameters,ax,...
      "Data",{ic_p},...
      "Meta",{struct('name',"Single Peak", ...
                     'discrete',false, ...
                     'color',[0,0,255]/255, ...
                     'thickness',0.01)}, ...
      "Legend",true,...
      "Regions",true,...
      "Region_labels",true,...
      "Title","Initial Conditions");

saveas(gcf,'single_peak_time_0.png')

disp('Visualized IC')

%% Simulate Continuous System (MVP)

[ps_t, ps_d] = cont_sudospec_etd_feuler(ic_p,parameters);

disp('Simulated MVP')

%% Visualize Ending Distribution

figure(2)
clf
ax = gca;
frame(x,parameters,ax,...
      "Data",{ps_d(end,:)},...
      "Meta",{struct('name',"Single Peak", ...
                     'discrete',false, ...
                     'color',[0,0,255]/255, ...
                     'thickness',0)}, ...
      "Legend",true,...
      "Regions",true,...
      "Title","Time=10");

saveas(gcf,'single_peak_end_time.png')

disp("Visualized Final Distribution")

%% View Animation

parameters.fr    = 50; 
parameters.pr    = 0.002;

figure(3)
clf

ax = gca;

animate(x,parameters,ax, ...
    Data={ps_d}, ...
    Meta={struct('name','Single Peak', ...
                'discrete',false, ...
                'color',[0,0,0], ...
                'thickness',parameters.eps)}, ...
    Time=ps_t, ...
    Title="Test", ...
    Legend=true, ...
    Regions=true, ...
    Region_labels=true);

disp("Finished Animating")

%% Monitor Broadening

% Helper function for comparing the state at a given time to a gaussian
function return_data = gaussian_compare(x,y,animate)

    % if desired, set up environment for animation
    if animate
        figure(4)
        clf
        hold on
    end

    % Find location of maximum value 
    [M, m] = max(y);

    % Compute center of domain
    d = round(length(y)/2);

    % Shift data to be centered on domain (otherwise periodic boundaries
    % make regression more difficult)
    if d > m
        y = [y(end-(abs(m-d)-1):end) y(1:end-abs(m-d))];
    end
    if m > d
        y = [y(abs(m-d)+1:end) y(1:abs(m-d))];
    end

    % Find region of positivity
    be = find(y>exp(-10),1,"first");
    en = find(y>exp(-10),1,"last");

    % Restrict to positive region and take logarithmic. In the case of a
    % gaussian, this should yield a quadratic
    f = y(be:en);
    lf = log(f);

    % Built-in Matlab polynomial fit
    [pcf, rrr] = polyfit(x(be:en),lf,2);

    % if desired, plot comparison of data to fit
    if animate
        plot(x(be:en),lf,LineWidth=4,DisplayName="Data")
        plot(x(be:en),polyval(pcf,x(be:en)),LineWidth=2,DisplayName="Fit")
        legend();
    end

    % return standard-deviation fitted gaussian, as well as r-squared value
    return_data = [1/sqrt(-pcf(1)), rrr.rsquared];

end

s2 = size(ps_d);
s2 = s2(1);

track_sd = zeros([1,s2]);
track_r2 = zeros([1,s2]);

% Study fitted gaussian at every timestep
for k = 1:s2
    out = gaussian_compare(x,ps_d(k,:),false);
    track_sd(k) = out(1);
    track_r2(k) = out(2);
    fprintf("Done %d out of %d\n",k,s2)
end

disp("Collected Spread Data")

%% Compare to Heat Propagation

% Square dilation. In case of gaussian under heat dynamics, this should
% yield a line.
track_sd_2 = track_sd.^2;

% Built-in Matlab polynomial fit
[pcf, rrr] = polyfit(ps_t,track_sd_2,1);
r2 = rrr.rsquared;

fprintf("Fit value of %f with r-squared of %f\n",pcf(1),r2);
fprintf("This suggests an epsilon of %f\n",pcf(1)/start_width);

disp("Fit Spread")

%% Visualize Spreading

figure(5)
clf
hold on

% Plot gaussian bandwidth over time
plot(ps_t,track_sd,LineWidth=4,DisplayName="Data")
%plot(ps_t,sqrt(pcf(1)*ps_t),LineWidth=2,DisplayName="Fit")
legend();
title("Broadening of Gaussian")

saveas(gcf,'broadening.png')

bad_r2 = min(track_r2);
fprintf("Worst regression was %f\n",bad_r2)

disp("Tracked Spread")