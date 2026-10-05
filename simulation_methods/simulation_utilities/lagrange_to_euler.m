function return_data = lagrange_to_euler(lag_coord, eul_coord, mass)
%LAGRANGE_TO_EULER converts lagrangian coordinates to the original eulerian
%coordinates in which initial data were given. This is purely for
%presentation, it does not impact the algorithm.
%
%last updated 09/23/26 by Adam Petrucci

    % Periodic extension (for easy vectorization)
    d_ext = [lag_coord(end)-1, lag_coord, lag_coord(1)+1];
    ic_ext = [eul_coord,1];
    mass_ext = [mass(end), mass];
    cum_mass = [0, cumsum(mass_ext)];

    % Compute density on each cell
    density_lag = mass_ext ./ diff(d_ext);

    % linear interpolation of mass in Eulerian cells, leveraging
    % monotonicity of vectors
    interp_mass = zeros(size(ic_ext));
    k = 1;
    for i = 1:length(ic_ext)
        while k < length(mass_ext) && ic_ext(i) >= d_ext(k+1)
              k = k + 1;
        end
        interp_mass(i) = cum_mass(k) ...
                + density_lag(k) * (ic_ext(i) - d_ext(k));
    end

    % density on Eulerian mesh
    return_data = diff(interp_mass) ./ diff(ic_ext);

end