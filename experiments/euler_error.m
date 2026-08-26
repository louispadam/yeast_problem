function return_data = euler_error(ic_dists,params,options)
%EULER_ERROR computes Wasserstein error in solutions of NODE simulated
% with an Euler scheme. This leverages the fact that the proof scheme
% computes exact (up to machine error) solutions. The method is designed to
% run across a vector of number of particles (N), many geometries
% (characterized by delta), a vector of stepsizes, and a cell array of
% distribution functions to sample initial conditions from. The geometry in
% question fixes |R|=|S|=0.3.
%
%last updated 07/31/26 by Adam Petrucci
arguments (Input)
    ic_dists              % initial distributions to sample from
    params                % geometric parameters
end
arguments (Input)
    options.X_data = linspace(0,1,2^11)  % spatial discretization (for
                                         % sampling initial conditions)
    options.N = [6000]                   % range of N to explore
    options.Delta = [0.2]                % range of |delta| to explore
    options.Trials = 10                  % number of trials to average
    options.Stepsizes = 0.002 * (1:1:50) % timestep for Euler scheme
    options.Runtime = 50                 % length of simulation time
    options.Figure = -1;                 % produce figure
                                         % Default: no figure
end
arguments (Output)
    return_data
end

    %****************************
    % Collect Inputs
    %****************************

    x = options.X_data;
    N = options.N;
    delta = options.Delta;
    trials = options.Trials;
    stepsizes = options.Stepsizes;
    runtime = options.Runtime;

    error = zeros([length(ic_dists),length(N),length(delta),length(stepsizes)]);

    %****************************
    % Run Experiment
    %****************************

    for ic_ind = 1:length(ic_dists)

        for N_ind = 1:length(N)

            for delta_ind = 1:length(delta)

                % set the geometry
                del = delta(delta_ind);
                params.s1 = del/2;
                params.s2 = del/2 + 0.3;
                params.r1 = 0.7 - del/2;
                params.r2 = 1 - del/2;

                error_temp = zeros([1,trials]);

                % Simulate this parameter set with both proof and euler
                % schemes as many times as 'trials' and store average error
                for dt_ind = 1:length(stepsizes)

                    for i = 1:trials

                        % Sample initial conditions
                        ic_p = sort(sample(x,ic_dists{ic_ind},N(N_ind)));

                        % Simulate with proof scheme
                        [time_new, data_new, ~] = particle_proof(ic_p,params,...
                                                     Update=false, ...
                                                     EndTime=runtime);
                        new_scheme = data_new(end,:);

                        % Simulate with Euler scheme
                        [time_old, data_old, ~] = particle_feuler(ic_p,params, ...
                                                    Update=false, ...
                                                    Timestep=stepsizes(dt_ind), ...
                                                    EndTime=runtime);
                        [~, ind] = min(abs(time_old - time_new(end)));
                        old_scheme = data_old(ind,:);

                        error_temp(i) = metric_wasserstein1(old_scheme,new_scheme);

                    end % of trials for this parameter set

                    % Compute average error
                    error(ic_ind,N_ind,delta_ind,dt_ind) = mean(error_temp);

                end % of stepsize loop

            end % of delta loop

        end % of N loop

    end % of ic loop

    return_data = error;

    %****************************
    % Generate Figure
    %****************************

    fig_num = options.Figure;

    if fig_num > 0

        error_frame = figure(fig_num);
        clf(error_frame);
        ax = axes(error_frame);

        plot(ax,stepsizes,squeeze(error(1,1,1,:)),...
                Color='black',...
                LineWidth=1.5)

        title(ax,'Error in Euler Scheme')
        xlabel(ax,'Stepsize')
        ylabel(ax,'Wasserstein Error')

        % Save figure
        timestamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
        filename = ['euler_error_' timestamp '.png'];
        saveas(error_frame,filename)

    end

end