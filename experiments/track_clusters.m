function return_data = track_clusters(time_vec,data,options)
%TRACK_CLUSTERS
%
%last updated 08/22/26 by Adam Petrucci
arguments (Input)
    time_vec       % initial conditions
    data           % parameters for simulation
end
arguments (Input)
    options.Modes = 7
    options.Discrete = true
end
arguments (Output)
    return_data
end

    to_track = options.Modes;

    track_estuff = zeros([length(time_vec),to_track]);
 
    for i = 1:length(time_vec)
        estuff = quantify_clustering(squeeze(data(i,:)),...
                                     "Discrete",options.Discrete,...
                                     "Modes",options.Modes);
        track_estuff(i,:) = abs(estuff);
    end

    return_data = track_estuff;

end