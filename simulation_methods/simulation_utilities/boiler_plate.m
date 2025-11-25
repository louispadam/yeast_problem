function return_data = boiler_plate()

    % Generate parameter structure
    parameters = struct(...
        's1',0, ...        % beginning of signalling region
        's2',0.3, ...      % end of signalling region
        'r1',0.4, ...      % beginning of responsive region
        'r2',0.5, ...      % end of responsive region
        'dt',0.001, ...    % simulation time step
        'tfin',100, ...    % ending time
        'del',0.1, ...     % sharpness of hyperbolic cutoffs
        'eps',0.025, ...   % diffusive coefficient
        'alph',0.9, ...    % linear inhibition term
        'fr',2, ...        % frame rate
        'pr',1, ...        % pause rate
        'ct',@(x) 1, ...   % cutoff function
        'm_sz',1, ...      % max vector size
        'update',false);   % print progress of simulation

    return_data = parameters;

end