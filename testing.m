%% testing.m
%
% This file includes a collection of tests for numerical methods. One
% should be able to roughly evaluate the qualitative accuracy of a
% numerical simulation of the yeast problem using these tests.
%
% last updated 11/13/25 by Adam Petrucci

% Generate parameter structure
parameters = boiler_plate();

% Set parameters
parameters.s1     = 0.1;
parameters.s2     = 0.4;
parameters.r1     = 0.6;
parameters.r2     = 0.7;
parameters.dt     = 0.001;
parameters.tfin   = 10;
parameters.del    = 0.05;
parameters.eps    = 0.05;
parameters.alph   = 0.9;
parameters.fr     = 50; 
parameters.pr     = 0.002;
parameters.ct     = ct_sharp(parameters.del);
parameters.m_sz   = 2^15;
parameters.update = true;

to_test = @(x,y) cont_sudospec_beuler(x,y);

% Test a narrow Gaussian
test_gaussian(to_test, parameters);

% Test the stationary solution
%test_stationary(to_test, parameters);

% Test smoothed-out stationary solution
%test_stationary_smooth(to_test, parameters);

%-HELPER-FUNCTIONS-----------------------------------------------------

function test_gaussian(to_test, parameters)

    % Setup domain
    n = 11;
    x = linspace(0,1,2^n);

    % Generate initial conditions
    ex = @(z) exp(-(z-0.5).^2/0.0001);
    ic = ex(x)/lp_integrate(x,ex(x),1);

    run_test(x,ic,to_test,parameters);

end

function test_stationary(to_test, parameters)

    % Setup domain
    n = 11;
    x = linspace(0,1,2^n);

    % Generate initial conditions
    ic = stationary_soln_v(x,parameters);

    run_test(x,ic,to_test,parameters);

end

function test_stationary_smooth(to_test,parameters)

    % Setup domain
    n = 11;
    x = linspace(0,1,2^n);

    % Generate initial conditions
    ic = stationary_soln_v(x,parameters);
    ex = @(z) exp(-(z-0.5).^2/0.0001);
    e = ex(x)/lp_integrate(x,ex(x),1);
    tempe = conv(ic,e,'same')/length(x);
    ic(100:end-100) = tempe(100:end-100);

    run_test(x,ic,to_test,parameters);

end

function run_test(x,ic,to_test,parameters)

    % Visualize Initial Conditions

    figure(1)
    clf
    ax = gca;
    frame(x,parameters,ax,...
        "Data",{ic},...
        "Meta",{struct('name',"Initial Distribution", ...
                       'discrete',false, ...
                       'color',[0,0,255]/255, ...
                       'thickness',0.01)}, ...
        "Legend",true,...
        "Regions",true,...
        "Region_labels",true,...
        "Title",["Time: " 0]);

    disp('Visualized IC')

    %---------------------------------------------------------------

    % Simulate Continuous System (MVP)

    [ps_t, ps_d] = to_test(ic,parameters);

    disp('Simulated MVP')

    %---------------------------------------------------------------

    % Visualize Ending Distribution

    figure(2)
    clf
    ax = gca;
    frame(x,parameters,ax,...
        "Data",{ps_d(end,:)},...
        "Meta",{struct('name',"Ending Distribution", ...
                       'discrete',false, ...
                       'color',[0,0,255]/255, ...
                       'thickness',0)}, ...
        "Legend",true,...
        "Regions",true,...
        "Title",["Time: " parameters.tfin]);

    disp("Visualized Final Distribution")

end