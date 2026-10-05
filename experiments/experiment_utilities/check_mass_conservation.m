function return_data = check_mass_conservation(x,time,data,options)
%CHECK_MASS_CONSERVATION checks the mass at every timestep of a simulation,
%and can plot the results
%
%last updated 09/24/25
arguments (Input)
    x         % discretization of domain
    time      % simulation timesteps
    data      % simulation data
end
arguments (Input)
    options.FigNum = 0  % whether or not to plot (default is not)
end
arguments (Output)
    return_data    % mass at every timestep
end

    % Prepare objects
    fignum = options.FigNum;
    track_mass = zeros(size(time));

    % Compute mass (on periodic domain) at every timestep
    for k = 1:length(times)
        dd = squeeze(data(k,:));
        track_mass(k) = 0.5 * sum([(dd(2:end) + dd(1:end-1)) .* (x(2:end) - x(1:end-1)),...
                                   (dd(1) + dd(end)) * mod(x(1) - x(end),1)]);
    end

    % Code for plotting, if desired
    if fignum

        % Set up figure
        mass_figure = figure(fignum);
        clf(mass_figure);
        ax = axes(mass_figure);

        plot(ax,1:length(time),track_mass)

        % Print basic statistics
        fprintf(['Max mass     : %.16f\n', ...
                 'Average mass : %.16f\n', ...
                 'Min mass     : %.16f\n'], ...
                max(track_mass), mean(track_mass), min(track_mass));

    end

    return_data = track_mass;

end