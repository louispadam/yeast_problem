function return_data = process_pseudospectral(domain,time,data)
%PROCESS_PSEUDOSPECTRAL takes data from the pseudospectral
%method and transforms it for the purposes of the animate or frame
%graphical methods
arguments
    domain (1,:)    % discretization of domain
    time (1,:)      % time steps in simulation
    data (:,:)      % results of a simulation via pseudospectral
end

    % Collect Inputs
    dom = domain;

    % Add dimension to data object
    ps_f = zeros([1,length(time),length(dom)]);
    for k = 1:length(time)
        ps_f(1,k,:) = data(k,:);
    end

    return_data = ps_f;
    
end