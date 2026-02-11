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

%% Generate Continuous IC

% Setup domain
n = 11;
x = linspace(0,1,2^n);

% Setup initial: a k-peak wave with some thickness at the base
k = 4;
ic_s = @(z) (1-cos(2*pi*k*z))+0.5;
s_n = lp_integrate(x,ic_s(x),1);
ic = @(z) ic_s(z)/s_n;

ic_c = ic(x);

disp("Generated Continuous IC")

%% Recover Particle IC

load("02_03_ic.mat","ic_p")

disp("Recovered Particle IC")

%% Recover Experiment Data

load("02_03_data_full.mat","exp_colec")

disp("Recovered Experiment Data")

%%

function return_data = test_peaks(x,d,options)
arguments (Input)
    x       % x data
    d       % y data
end
arguments (Input)
    options.mass_tol double = 0.2
    options.height_tol double = 1
end

    return_data = options.mass_tol > lp_integrate(x,min(d,options.height_tol*ones(1,length(d))),1);

end

%%

found_clusters = zeros(4,0);
total_count = 0;

for i1 = 1:length(exp_colec)
    for i2 = 1:length(exp_colec{i1})
        for i3 = 1:length(exp_colec{i1}{i2})
            for i4 = 1:length(exp_colec{i1}{i2}{i3})
                if test_peaks(x,exp_colec{i1}{i2}{i3}{i4}{2}{2})
                    found_clusters(:,end+1) = [i1;i2;i3;i4];
                end
                total_count = total_count + 1;
            end
        end
    end
end

fprintf("Found %d candidates of %d experiments\n",size(found_clusters,2),total_count);

%%

function [return_bool,return_vals] = compare_peaks(x,d1,d2,th,options)
arguments (Input)
    x
    d1
    d2
    th
end
arguments (Input)
    options.mph double = 1
    options.mpd double = 0
    options.ha logical = false
end

    if length(x) ~= length(d1)
        d1 = fatten_points_polynomial(x,d1,th);
    end

    if length(x) ~= length(d2)
        d2 = fatten_points_polynomial(x,d2,th);
    end

    peaks1 = findpeaks(d1,x,"MinPeakHeight",options.mph,...
                            "MinPeakDistance",options.mpd);
    peaks2 = findpeaks(d2,x,"MinPeakHeight",options.mph,...
                            "MinPeakDistance",options.mpd);

    if options.ha
        peaks1(peaks1 < max(peaks1)/2) = [];
        peaks2(peaks2 < max(peaks2)/2) = [];
    end

    return_bool = (length(peaks1) == length(peaks2));
    return_vals = [length(peaks1),length(peaks2)];

end

%%

peak_differences = zeros(4,0);

for i = 1:size(found_clusters,2)

    p_i = found_clusters(:,i);

    c_test = exp_colec{p_i(1)}{p_i(2)}{p_i(3)}{p_i(4)}{2}{2};
    p_test = exp_colec{p_i(1)}{p_i(2)}{p_i(3)}{p_i(4)}{2}{1};

    check = compare_peaks(x,c_test,p_test,parameters.del,"mph",1,"mpd",0.1,"ha",true);

    if ~check(1)
        peak_differences(:,end+1) = p_i;
    end

end

fprintf("Found %d candidates of %d experiments\n",size(peak_differences,2),size(found_clusters,2));

%%

p_i = peak_differences(:,20);

c_test = exp_colec{p_i(1)}{p_i(2)}{p_i(3)}{p_i(4)}{2}{2};
p_test = exp_colec{p_i(1)}{p_i(2)}{p_i(3)}{p_i(4)}{2}{1};

[c_bool,c_vals] = compare_peaks(x,c_test,p_test,parameters.del,"mph",1,"mpd",0.1);
c_vals

%%

check = 2;
p_i = peak_differences(:,check);

p_d = exp_colec{p_i(1)}{p_i(2)}{p_i(3)}{p_i(4)}{1};
parameters.s1 = p_d(1);
parameters.s2 = p_d(2);
parameters.r1 = p_d(3);
parameters.r2 = p_d(4);

c_test = exp_colec{p_i(1)}{p_i(2)}{p_i(3)}{p_i(4)}{2}{2};
p_test = exp_colec{p_i(1)}{p_i(2)}{p_i(3)}{p_i(4)}{2}{1};

test_frame = figure(2);
clf(test_frame);
ax = axes(test_frame);

frame(x,parameters,ax,...
      "Data",{c_test,p_test},...
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
      "Title",sprintf("Testing Candidate %d",check));

disp("Visualized Test")

%%

check = 2;
p_i = peak_differences(:,check);

p_d = exp_colec{p_i(1)}{p_i(2)}{p_i(3)}{p_i(4)}{1};
parameters.s1 = p_d(1);
parameters.s2 = p_d(2);
parameters.r1 = p_d(3);
parameters.r2 = p_d(4);

%%

num_particles = 1000:5:1200;

pks = zeros(1,length(num_particles));

for i = 1:length(num_particles)

    ic_p = sort(sample(x,ic,num_particles(i)));
    [t_p, d_p] = particle_feuler(ic_p,parameters);
    f_p = d_p(end,:);

    [c_bool,c_vals] = compare_peaks(x,exp_colec{p_i(1)}{p_i(2)}{p_i(3)}{p_i(4)}{2}{2},...
                                    f_p,parameters.del,"mph",1,"mpd",0.1);
    pks(i) = c_vals(1)/c_vals(2);

    fprintf("Done %d of %d\n",i,length(pks));

end

%%

pro_frame = figure(2);
clf(pro_frame);
ax = axes(pro_frame);

plot(ax,num_particles,pks)
ylim(ax,[0,3])

%%

one_frame = figure(2);
clf(one_frame);
ax = axes(one_frame);

frame(x,parameters,ax,...
      "Data",{exp_colec{p_i(1)}{p_i(2)}{p_i(3)}{p_i(4)}{2}{2},f_p},...
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
      "Title",sprintf("%d particles",num_particles));

%%

1:2:10