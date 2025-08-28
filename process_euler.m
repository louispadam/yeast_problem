function return_data = process_euler(domain,time,data,std_dev)
%PROCESS_EULER takes data from the forward_euler or forward_euler_noise
%methods and transforms it for the purposes of the animate or frame
%graphical methods
arguments
    domain (1,:)    % discretization of domain
    time (1,:)      % time steps in simulation
    data (:,:)      % results of a simulation via euler
    std_dev         % thickening parameter
end

    % Collect Inputs
    dom = domain;
    std = std_dev;
    
    % Process NODE data by thickening points, and add dimension to data
    % object
    fe_f = zeros([1,length(time),length(dom)]);
    for k = 1:length(time)
        fe_f(1,k,:) = fatten_points_polynomial(dom,data(k,:),std);
    end

    return_data = fe_f;
    
end