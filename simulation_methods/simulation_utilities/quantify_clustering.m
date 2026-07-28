function return_data = quantify_clustering(state,options)
arguments (Input)
    state
end
arguments (Input)
    options.Modes = 5;
end

    Modes = options.Modes;

    % convert to radians
    angles = 2*pi*state;

    ft = exp(1i*angles);
    ftk = ft;

    modes_vec = ones([1,Modes]);

    for k = 1:Modes
        modes_vec(k) = mean(ftk);
        ftk = ftk .* ft;
    end

    return_data = modes_vec;

end