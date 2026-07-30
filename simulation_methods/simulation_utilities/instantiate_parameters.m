function return_data = instantiate_parameters()

    % Generate parameter structure
    parameters = struct(...
        's1',0, ...        % beginning of signalling region
        's2',0.3, ...      % end of signalling region
        'r1',0.4, ...      % beginning of responsive region
        'r2',0.5, ...      % end of responsive region
        'eps',0.025, ...   % diffusive coefficient
        'alph',0.9, ...    % linear inhibition term
        'ct',@(x) 1);      % cutoff function

    return_data = parameters;

end