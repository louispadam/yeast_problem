function return_data = track_clusters(time_vec,data,options)
%TRACK_CLUSTERS estimates the number of clusters through a simulation of
%the yeast system (particle or continuum).
%
%last updated 11/13/26 by Adam Petrucci
arguments (Input)
    time_vec       % initial conditions
    data           % parameters for simulation
end
arguments (Input)
    options.Modes = 7
    options.Type = "Discrete"
    options.IC_Mass = []
end
arguments (Output)
    return_data
end

    to_track = options.Modes;

    if strcmp(options.Type,"Lagrangian") 
        if isempty(options.IC_Mass)
            error("Specified Lagrangian format but did not provide " + ...
                  "mass vector.")
        end
        ic_mass = options.IC_Mass;
        to_send = @(z) process_lagrangian(z,ic_mass);

    elseif strcmp(options.Type,"Eulerian")
        mesh_size = size(data,2);
        mesh = linspace(0,1-1/mesh_size,mesh_size);
        to_send = @(z) [mesh;z];

    elseif strcmp(options.Type,"Discrete")
        num_par = size(data,2);
        mass = ones([1,num_par])/num_par;
        to_send = @(z) [z;mass];

    else
        error("Did not recieve acceptable data format.")
    end

    track_estuff = zeros([length(time_vec),to_track]);
 
    for i = 1:length(time_vec)
        dd = squeeze(data(i,:));
        estuff = quantify_clustering(to_send(dd),...
                                     "Modes",options.Modes);
        track_estuff(i,:) = abs(estuff);
    end

    return_data = track_estuff;

end

function return_data = process_lagrangian(dd,ic_mass)
    mesh = dd(2:end);
    mass = circshift(ic_mass,-dd(1));
    return_data = [mesh;mass];
end