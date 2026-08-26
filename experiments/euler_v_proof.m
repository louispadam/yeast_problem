function [return_proof_clocks, return_euler_clocks, return_error] = ...
         euler_v_proof(ic_dists, params, options)
%EULER_V_PROOF compares runtimes of NODE simulations with an Euler scheme
% and the proof scheme. The method is designed to
% run across a vector of number of particles (N), many geometries
% (characterized by delta), and a cell array of distribution functions to
% sample initial conditions from. The geometry in question fixes
% |R|=|S|=0.3.
%
%last updated 07/31/26 by Adam Petrucci
arguments (Input)
    ic_dists                        % initial distributions to sample from
    params                          % geometric parameters
end
arguments (Input)
    options.X_data = linspace(0,1,2^11)  % spatial discretization (for
                                         % sampling initial conditions)
    options.N = 2000 * (1:1:10)          % range of N to explore
    options.Delta = 0.05 * (1:1:8)       % range of |delta| to explore
    options.Trials = 10                  % number of trials to average
    options.Stepsize = 0.01              % timestep for Euler scheme
    options.Runtime = 50                 % length of simulation time
    options.Figure = -1;                 % produce figure
                                         % Default: no figure
end
arguments (Output)
    return_proof_clocks
    return_euler_clocks
    return_error
end

    %****************************
    % Collect Inputs
    %****************************

    x = options.X_data;
    N = options.N;
    delta = options.Delta;
    trials = options.Trials;
    stepsize = options.Stepsize;
    runtime = options.Runtime;

    % Instantiate storage objects
    proof_clocks = zeros([length(ic_dists),length(N),length(delta)]);
    euler_clocks = zeros([length(ic_dists),length(N),length(delta)]);
    error = zeros([length(ic_dists),length(N),length(delta)]);

    %****************************
    % Run Experiment
    %****************************

    for ic_ind = 1:length(ic_dists)

        for N_ind = 1:length(N)

            for delta_ind = 1:length(delta)

                % Set geometry
                del = delta(delta_ind);
                params.s1 = del/2;
                params.s2 = del/2 + 0.3;
                params.r1 = 0.7 - del/2;
                params.r2 = 1 - del/2;

                proof_clocks_temp = zeros([1,trials]);
                euler_clocks_temp = zeros([1,trials]);
                error_temp = zeros([1,trials]);

                % Simulate this parameter set with both proof and euler
                % schemes as many times as 'trials' and store averages
                for i = 1:trials

                    % Sample initial conditions
                    ic_p = sort(sample(x,ic_dists{ic_ind},N(N_ind)));

                    % Simulate with proof scheme
                    [time_new, data_new, clock_new] = particle_proof(ic_p,params,...
                                                          Update=false, ...
                                                          EndTime=runtime);
                    new_scheme = data_new(end,:);

                    % Simulate with Euler scheme
                    [time_old, data_old, clock_old] = particle_feuler(ic_p,params, ...
                                                    Update=false, ...
                                                    Timestep=stepsize, ...
                                                    EndTime=runtime);
                    [~, ind] = min(abs(time_old - time_new(end)));
                    old_scheme = data_old(ind,:);

                    proof_clocks_temp(i) = clock_new;
                    euler_clocks_temp(i) = clock_old;

                    error_temp(i) = metric_wasserstein1(old_scheme,new_scheme);
                
                end % of trials for this parameter set

                % Compute average error
                proof_clocks(ic_ind,N_ind,delta_ind) = mean(proof_clocks_temp);
                euler_clocks(ic_ind,N_ind,delta_ind) = mean(euler_clocks_temp);
                error(ic_ind,N_ind,delta_ind) = mean(error_temp);

            end % of delta loop
    
        end % of N loop

    end % of ic loop

    return_proof_clocks = proof_clocks;
    return_euler_clocks = euler_clocks;
    return_error = error;

    %****************************
    % Generate Figure
    %****************************

    fig_num = options.Figure;

    if fig_num > 0

        scheme_comparison = figure(fig_num);
        clf(scheme_comparison);
        ax = axes(scheme_comparison);

        c = turbo(32);

        blueColors = c(1:8,:);            % blue/cyan region
        redColors  = flip(c(25:32,:));    % red/orange region

        hold on

        for i = 1:length(delta)

            plot(ax,N,log(proof_clocks(1,:,1+mod(i,length(delta)))),...
                      Color=redColors(i,:),...
                      LineWidth=1.5,...
                      DisplayName=sprintf('\\delta=%.2f',delta(i)));

        end

        for i = 1:length(delta)
    
            plot(ax,N,log(euler_clocks(1,:,1+mod(i,length(delta)))),...
                      Color=blueColors(i,:),...
                      LineWidth = 1.5,...
                      DisplayName=sprintf('\\delta=%.2f',delta(i)));
        end

        hold off

        title('Simulation Times with Uniform Sampling');
        xlabel('Number of Particles')
        ylabel('Total Simulation Time (log scale)')

        % The legend for this plot is cumbersome, so what follows is a
        % manual one

        % Create key axes
        axKey = axes('Position',[0.6 0.15 0.30 0.15]);
        hold(axKey,'on')
        axis(axKey,'off')

        % Coordinates
        xVals = 1:8;
        yBlue = 2;
        yRed  = 1;

        % Label for the row-header column
        text(axKey,0.2,3.1,sprintf('\\bf\\delta\\rm x 10^{-2} = '), ...
            'HorizontalAlignment','right');

        % Column headers (gradient levels)
        for i = 1:8
            text(axKey,xVals(i),3.0,sprintf('%.0f',100*delta(i)), ...
                'HorizontalAlignment','center', ...
                'FontWeight','bold');
        end

        % Row labels (classes)
        text(axKey,0.2,yBlue,'Euler Scheme', ...
            'HorizontalAlignment','right', ...
            'FontWeight','bold');

        text(axKey,0.2,yRed,'Proof Scheme', ...
            'HorizontalAlignment','right', ...
            'FontWeight','bold');

        % Draw blue entries
        for i = 1:8
            plot(axKey,[xVals(i)-0.25 xVals(i)+0.25], ...
                [yBlue yBlue], ...
                'Color',blueColors(i,:), ...
                'LineWidth',2);
        end

        % Draw red entries
        for i = 1:8
            plot(axKey,[xVals(i)-0.25 xVals(i)+0.25], ...
                [yRed yRed], ...
                'Color',redColors(i,:), ...
                'LineWidth',2);
        end

        xlim(axKey,[0 9]);
        ylim(axKey,[0.5 3.7]);

        % Save figure
        timestamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
        filename = ['euler_proof_comparison_' timestamp '.png'];
        saveas(error_frame,filename)

    end

end