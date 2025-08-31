function return_data = construct_em(x_data,y_data,num_points)
%CONSTRUCT_EM Sample a given probability distribution to construc an
%approximating empirical measure.
%
%last updated 08/30/25 by Adam Petrucci
arguments (Input)
    x_data (1,:) double % discretization of domain for y_data
    y_data (1,:) double % distribution to be sampled from
    num_points          % number of points in sample
end

    % Collect required inputs (for notational convenience)
    n = num_points;

    %****************************
    % Construct CDF of y_data
    %****************************
    dist_y = zeros([1,length(x_data)]);
    for k = 2:length(y_data)
        dist_y(k) = lp_integrate(x_data(1:k),y_data(1:k),1);
    end
    dist_y(1) = 0;
    dist_y(end) = 1;

    %****************************
    % Construct Inverse CDF of y_data
    %****************************
    inv_dist_y = zeros([1,length(x_data)]);
    inv_x = linspace(0,1,length(x_data));   % same discretization as x_data
    for k = 1:length(inv_x)
        inv_dist_y(k) = x_data(find(dist_y>=inv_x(k),1));
    end

    %****************************
    % Sample
    %****************************
    particles_uniform = rand([1,n]);
    particles_dist = zeros([1,n]);
    for k = 1:n % I wonder if this is trivially vectorizeable
        particles_dist(k) = inv_dist_y(find(inv_x>=particles_uniform(k),1));
    end

    %****************************
    % Separate Any Doubles (due to discretization) 
    %****************************
    for k = 1:length(particles_dist)

        % Find repeats
        repeats = find(particles_dist == particles_dist(k));

        % Fix repeats by spreading them out between within discretization
        if length(repeats) > 1
            current = particles_dist(k);
            replace = inv_x(2)-inv_x(1)*linspace(0,1,length(repeats)+2)+current;
            particles_dist(repeats) = replace(2:end-1);
        end
    end

    return_data = particles_dist;
end