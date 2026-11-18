function return_data = quantify_clustering(state,options)
%QUANTIFY_CLUSTERING effectively computes magnitude of the first few
%Fourier modes of an empirical measure. When the data is actually
%continuous, it treats it as equally-spaced point-masses. This is intended
%to be used as a proxy for cluster detection.
%
%last updated 11/13/26 by Adam Petrucci
arguments (Input) 
    state                     % state vector on which to detect clustering
end
arguments (Input)
    options.Modes = 5         % modes to measure
end

    Modes = options.Modes;

    mesh = state(1,:);
    mass = state(2,:);

    angles = 2*pi*mesh;
    F_angles = exp(1i*angles);
    F_kangles = F_angles;

    modes_vec = ones([1,Modes]);

    for k = 1:Modes
        modes_vec(k) = mean(mass .* F_kangles);
        F_kangles = F_kangles .* F_angles;
    end

    return_data = modes_vec;

end