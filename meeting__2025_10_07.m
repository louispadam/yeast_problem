%% 10/07/2025
%
% This script is for a research meeting with Keith Promislow. It is
% intended to compare simulation results now that the cutoff function is
% consistent across equations.
%
% last updated 10/06/25 by Adam Petrucci

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
parameters.ct    = ct_tanh_k(parameters.del);

disp('Set Parameters')

%% Generate Initial Conditions

% Setup domain
n = 15;
x = linspace(0,1,2^n);

% Setup continuous IC: a k-peak wave with some thickness at the base
k = 4;
ic_s = @(z) (1-cos(2*pi*k*z))+0.2;
s_n = lp_integrate(x,ic_s(x),1);
ic_d = @(z) ic_s(z)/s_n;

% Setup particle IC (by drawing from continuous IC)
trials = 3;
max_p = 1500;
ic_p = zeros([trials,max_p]);
for c = 1:trials
    ic_p(c,:) = sample(x,ic_d,max_p);
end

disp('Generated IC')

%% Simulate Continuous System (MVP)

[ps_t, ps_d] = pseudospectral(ic_d(x),parameters);

disp('Simulated MVP')

%% Process MVY Data

t_num = 6;
times = round(linspace(0,parameters.tfin,t_num));

to_save = ones([t_num,length(x)]);
for j = 1:length(times)
    t_j = find(ps_t>=times(j),1,"first");
    if isempty(t_j)
        t_j = length(ps_t);
    end
    to_save(j,:) = ps_d(t_j,:);
end
data_s = {to_save};

disp('Processed MVY data')

%% Simulate Particle System (NODE)

grain = 3;
p_num = round(linspace(max_p/grain,max_p,grain));

data_p = cell(grain,trials);

for m = 1:grain

    for c = 1:trials

        p = ic_p(c,1:p_num(m));
        [en_t, en_p] = forward_euler_noise(p,parameters);

        to_save = ones([t_num,p_num(m)]);
        for j = 1:t_num
            t_j = find(en_t>=times(j),1,"first");
            if isempty(t_j)
                t_j = length(ps_t);
            end
            to_save(j,:) = en_p(t_j,:);
        end

        data_p{m,c} = to_save; 

    end

    disp(['Done ',num2str(p_num(m)),'-particle trials']);

end

disp('Simulated NODE')

%% Graph Data

figure(1)
clf

tl_graph = tiledlayout(grain,t_num);

axs = gobjects(grain,t_num);

for np = 1:grain

    for grab_from = 1:t_num

        meta_d = cell(1,1+grain);

        meta_d{1} = struct('name',"Continuous", ...
            'color',[0,0,255]/255, ...
            'thickness',0);

        to_send = cell(1,1+trials);
        ar = data_s{1};
        to_send{1} = ar(grab_from,:);

        for i = 1:trials
            ar = data_p{np,i};
            to_send{1+i} = ar(grab_from,:);
            meta_d{1+i} = struct('name',sprintf('Trial %d',1+i), ...
                'color',[-1,-1,-1], ...
                'thickness',parameters.eps);
        end

        axs(np,grab_from) = nexttile;

        frame(x,parameters,axs(np,grab_from), ...
            Data = to_send, ...
            Meta = meta_d);

        axs(np,grab_from).XTickLabels = [];
        axs(np,grab_from).YTickLabels = [];
        xlabel(axs(np,grab_from), '');
        ylabel(axs(np,grab_from), '');

    end

end

sgtitle(figure(1),"Experiment",'FontWeight', 'bold');
%linkaxes(axes,'y')
xlabel(tl_graph,'Position');
ylabel(tl_graph,'Density')

% ====== Column Labels (Top) ======
labelWidth = 0.05;
labelHeight = 0.05;

for j = 1:length(times)
    pos = axs(end,j).Position;
    xpos = pos(1) + pos(3)/2 - labelWidth/2;
    ypos = 0.25*pos(2);
    annotation('textbox', ...
        [xpos, ypos, labelWidth, labelHeight], ...
        'String', times(j), ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold', ...
        'FontSize', 11);
end

% ====== Row Labels (Left) ======

for i = 1:grain
    pos = axs(i,1).Position;
    xpos = 0.8*pos(1);
    ypos = pos(2) + pos(4)/2 - 0.9*labelHeight;
    annotation('textbox', ...
        [xpos, ypos, labelWidth, labelHeight], ...
        'String', p_num(i), ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'Rotation', 90, ...
        'FontWeight', 'bold', ...
        'FontSize', 11);
end

annotation('textbox', ...
    [0.9*axs(end,1).Position(1),0.25*axs(end,1).Position(2), labelWidth, labelHeight], ...
    'String', 'Time', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'FontWeight', 'bold', ...
    'FontSize', 11);
annotation('textbox', ...
    [0.8*axs(1,1).Position(1),0.1*axs(1,1).Position(2), labelWidth, labelHeight], ...
    'String', 'Particles', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'Rotation', 90, ...
    'FontWeight', 'bold', ...
    'FontSize', 11);

disp('Graphed Initial Conditions')

%%

meta_d = cell(1,2);

meta_d{1} = struct('name',"Continuous", ...
    'color',[0,0,255]/255, ...
    'thickness',0);

to_send = cell(1,2);
ar = data_s{1};
to_send{1} = ar(1,:);

ar = data_p{np,1};
to_send{2} = ar(grab_from,:);
meta_d{2} = struct('name',sprintf('Trial %d',2), ...
    'color',[-1,-1,-1], ...
    'thickness',parameters.eps);

figure(3)
clf

aa = gca;

frame(x,parameters,aa, ...
    Data = to_send, ...
    Meta = meta_d);