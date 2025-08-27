function return_data = construct_em(x_data,y_data,n)

    dist_y = zeros([1,length(x_data)]);

    for k = 2:length(y_data)
        dist_y(k) = lp_integrate(x_data(1:k),y_data(1:k),1);
    end
    dist_y(1) = dist_y(2);
    dist_y(end) = 1;

    inv_dist_y = zeros([1,length(x_data)]);
    inv_x = linspace(0,1,length(x_data));

    for k = 1:length(inv_x)
        inv_dist_y(k) = x_data(find(dist_y>=inv_x(k),1));
    end

    particles_uniform = rand([1,n]);
    particles_dist_nu = zeros([1,n]);

    for k = 1:n % I wonder if this is trivially vectorized
        particles_dist_nu(k) = inv_dist_y(find(inv_x>=particles_uniform(k),1));
    end

    for k = 1:length(particles_dist_nu)
        repeats = find(particles_dist_nu == particles_dist_nu(k));
        if length(repeats) > 1
            current = particles_dist_nu(k);
            next_big = min(particles_dist_nu(particles_dist_nu-current>0));
            if isempty(next_big)
                next_big = 1;
            end
            next_small = max(particles_dist_nu(current-particles_dist_nu>0));
            if isempty(next_small)
                next_small = 0;
            end
            replace = (next_big-next_small)*linspace(0,1,length(repeats)+2)+next_small;
            particles_dist_nu(repeats) = replace(2:end-1);
        end
    end

    return_data = particles_dist_nu;

end