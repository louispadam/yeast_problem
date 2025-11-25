%% 11/26/2025
%
% This script is for a research meeting with Keith Promislow.
%
% last updated 11/24/25 by Adam Petrucci

%% Boilerplate

addpath(genpath(pwd));
parameters = boiler_plate();

disp('Ran Boilerplate')

%% Set up problem

parameters.s1     = 0.45;      % very small and separated regions should 
parameters.s2     = 0.5;       % have little impact on mechanics, in
parameters.r1     = 0.1;       % particular diffusion
parameters.r2     = 0.15;
parameters.dt     = 0.001;
parameters.tfin   = 10;
parameters.del    = 0.01;
parameters.eps    = 0.00;      % NO DIFFUSION
parameters.alph   = 0.00;      % NO INTERACTION
parameters.fr     = 1; 
parameters.pr     = 0.002;
parameters.ct     = ct_tanh_a(parameters.del);
parameters.m_sz   = 2^15;
parameters.update = true;

disp('Set Parameters')

%% Generate Distributions to Compare

% Setup domain
n = 11;
x = linspace(0,1,2^n);

% Set up initial distrubution: a narrow gaussian
start_width = 0.05;
ic_s = @(z) exp(-(z-0.5).^2/(start_width^2));
lp = lp_integrate(x,ic_s(x),1);

d1 = @(z) ic_s(z)/lp;
d2 = @(z) ic_s(z-0.4)/lp;

d1 = d1(x);
d2 = d2(x);

disp("Generated Distributions")

%%

disp('start')

hh = metric_wasserstein_new2(x,d1,d2);

hh

disp('test')

%% Visualize Initial Conditions

figure(1)
clf
ax = gca;
frame(x,parameters,ax,...
      "Data",{d1,d2},...
      "Meta",{struct('name',"Distribution 1", ...
                     'discrete',false, ...
                     'color',[0,0,255]/255, ...
                     'thickness',0.01) ...
              struct('name',"Distribution 2", ...
                     'discrete',false, ....
                     'color',[100,100,255]/255, ...
                     'thickness',0.01)}, ...
      "Legend",true,...
      "Regions",false,...
      "Region_labels",false,...
      "Title","Initial Conditions");

saveas(gcf,'single_peak_time_0.png')

disp('Visualized IC')

%% Comparing Distributions

was = metric_wasserstein_new(x,d1,x,d2,parameters,1);

fprintf('\nFinal distance of %f\n',was);

%%

X = x;
Y = x;
CX = d1;
CY = d2;

Z = union(X,Y);

 = stretch_values(x,z)
y_new = stretch_values(y,z)

function return_data = total_cost(z,x,y)
    total = 0;
    for k = 2:length(z)
        total = total + unit_cost(x(k),y(k))*(z(k)-z(k-1));
    end
    total = unit_cost(x(1),y(1))*(z(1)-z(end)+1);
    return_data = total;
end

function return_data = unit_cost(x,y,lambda)
    return_data =(abs(x-y)).^lambda;
end

function return_data = stretch_values(x,z)

    nz = length(z);

    zx = nan([nz,1]);
    zx(ismember(z,x)) = x;
    if isnan(zx(1))
        zx(1) = zx(find(~isnan(zx),1,"first"));
    end

    z1x = ones([nz,nz]) * diag(zx)';
    z1x(triu(true(size(z1x)),1)) = NaN;

    z2x = z1x .* triu(ones([nz, nz]))';

    z3x = max(z2x, [], 2)';

    return_data = z3x;

end

disp('run')

