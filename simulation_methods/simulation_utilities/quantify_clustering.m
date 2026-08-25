function return_data = quantify_clustering(state,options)
%QUANTIFY_CLUSTERING effectively computes magnitude of the first few
%Fourier modes of an empirical measure. When the data is actually
%continuous, it treats it as equally-spaced point-masses. This is intended
%to be used as a proxy for cluster detection.
%
%last updated 08/25/26 by Adam Petrucci
arguments (Input) 
    state                     % state vector on which to detect clustering
end
arguments (Input)
    options.Modes = 5         % modes to measure
    options.Discrete = true   % data is either discrete (particle) or continuous
end

    Modes = options.Modes;

    % convert to radians
    if options.Discrete
        angles = 2*pi*state;
        ft = exp(1i*angles);
        state = ones(size(state));
    else
        x = 2*pi*linspace(0,1-1/length(state),length(state));
        ft = exp(1i*x);
    end

    ftk = ft;

    modes_vec = ones([1,Modes]);

    for k = 1:Modes
        modes_vec(k) = mean(state .* ftk);
        ftk = ftk .* ft;
    end

    return_data = modes_vec;

end