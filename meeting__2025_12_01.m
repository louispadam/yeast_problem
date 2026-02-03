%% 12/01/2025
%
% This script is for a research meeting with Keith Promislow. It is
% intended to exhibit a method of measuring diffusion as well as the issues
% with our diffusion code.
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
parameters.tfin   = 10;
parameters.del    = 0.01;
parameters.eps    = 0.005;
parameters.alph   = 0.0;
parameters.fr     = 20; 
parameters.pr     = 0.002;
parameters.ct     = ct_tanh_k(parameters.del);
parameters.m_sz   = 2^15;
parameters.update = true;

disp('Set Parameters')

%% Generate Initial Conditions

% Setup domain
n = 15;
x = linspace(0,1,2^n);

% Set up initial distrubution: a narrow gaussian
start_width = 0.005;
ic_s = @(z) exp(-z.^2/(start_width^2));
lp = lp_integrate(x,ic_s(x),1);
ic = @(z) ic_s(z)/lp;

ic_p = ic(x-0.5);

disp("Generated Initial Conditions")

%% Visualize Initial Conditions

initial_frame = figure(1);
clf(initial_frame);
ax = axes(initial_frame);

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

disp('Visualized IC')

%% Simulate Continuous System (MVP)

[t_c, d_c] = cont_sudospec_etd_feuler(ic_p,parameters);

disp('Simulated MVP')

%% Visualize Ending Distribution

final_frame = figure(2);
clf(final_frame);
ax = axes(final_frame);

frame(x,parameters,ax,...
      "Data",{d_c(end,:)},...
      "Meta",{struct('name',"Single Peak", ...
                     'discrete',false, ...
                     'color',[0,0,255]/255, ...
                     'thickness',0)}, ...
      "Legend",true,...
      "Regions",true,...
      "Title","Time=10");

disp("Visualized Final Distribution")

%% View Animation

animation_fig = figure(3);
clf(animation_fig);
ax = axes(animation_fig);

animate(x,parameters,ax, ...
    Data={d_c}, ...
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
function return_data = gaussian_compare(x,y)

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

    % Restrict to positive region and take logarithm. In the case of a
    % gaussian, this should yield a quadratic
    f = y(be:en);
    lf = log(f);

    figure(10)
    clf
    plot(x(be:en)-0.5,lf)

    % Built-in Matlab polynomial fit
    %[pcf, rrr] = polyfit(x(be:en),lf,2);

    return_data = polyfit(x(be:en)-0.5,lf,2);

    % return standard-deviation fitted gaussian, as well as r-squared value
    %return_data = [1/sqrt(-pcf(1)), rrr.rsquared];

end

% The above function is good, but I think we can be more efficient if we
% already expect the Gaussian to relate to initial conditions
function return_data = ic_compare(x,y,ic_sigma,ic_coeff,t)

    % Find location of maximum value 
    [M, m] = max(y);

    % Compute center of domain
    d = round(length(y)/2);

    coeff = M/ic_coeff;
    beta = 2*pi/coeff;
    eps = (beta^2-1)/(2*ic_sigma*t);

    return_data = eps;

end

%iterates = length(d_t);
iterates = 180;

track_eps = zeros([1,iterates]);

ic_max = max(ic_p);
ic_sigma = 1/start_width^2;

% Study fitted gaussian at every timestep
for k = 1:180
    out = ic_compare(x,d_c(k,:),ic_sigma,ic_max,t_c(k));
    track_eps(k) = out;
    %fprintf("Done %d out of %d\n",k,iterates)
end

out = gaussian_compare(x,d_c(2,:))

disp("Collected Spread Data")

%%

figure(10)
clf
hold on

%eps = track_eps(2)
%beta = sqrt(4*ic_sigma*eps*t_c(2)+1)
%comp = ic_max*2*pi/eps*ic(x/beta-0.5)+1;

%plot(x,comp,'k');
%plot(x,d_c(2,:),'b')

plot(100:180,track_eps(100:180))

%% Compare to Heat Propagation

% Square dilation. In case of gaussian under heat dynamics, this should
% yield a line.
track_sd_2 = track_beta.^2;

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
plot(ps_t,track_beta,LineWidth=4,DisplayName="Data")
plot(ps_t,sqrt(pcf(1)*ps_t),LineWidth=2,DisplayName="Fit")
legend();
title("Broadening of Gaussian")

saveas(gcf,'broadening.png')

bad_r2 = min(track_norm);
fprintf("Worst regression was %f\n",bad_r2)

disp("Tracked Spread")