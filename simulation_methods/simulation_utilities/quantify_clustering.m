function return_data = quantify_clustering(state,options)
arguments (Input)
    state
end
arguments (Input)
    options.Modes = 5
    options.Discrete = true
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