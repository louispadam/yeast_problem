function return_data = lagrange_to_euler(lag_coord, eul_coord, mass, options)
%LAGRANGE_TO_EULER converts lagrangian coordinates to eulerian coordinates.
%
%last updated 10/08/26 by Adam Petrucci
arguments (Input)
    lag_coord
    eul_coord
    mass
end
arguments (Input)
    options.Rho = []
    options.Atom = []
end
arguments (Output)
    return_data
end

    rho = options.Rho;
    atom = options.Atom;
    if isempty(rho)
        [lag_coord, rho, atom] = split_sing_cont(lag_coord,mass);
    elseif isempty(atom)
        atom = zeros(size(rho));
    end

    % Continuous cell masses
    width = [diff(lag_coord), lag_coord(1)+1-lag_coord(end)];
    cont_mass = rho .* width;

    % Periodic extensions
    d_ext    = [lag_coord(end)-1, lag_coord, lag_coord(1)+1];
    rho_ext  = [rho(end), rho, rho(1)];
    mass_ext = [cont_mass(end), cont_mass, cont_mass(1)];
    atom_ext = [atom(end), atom, atom(1)];
    ic_ext = [eul_coord,1];

    cum_mass = [0,cumsum(mass_ext + atom_ext)];

    % linear interpolation of mass in Eulerian cells, leveraging
    % monotonicity of vectors
    interp_mass = zeros(size(ic_ext));
    k = 1;
    for i = 1:length(ic_ext)
        while k < length(d_ext) && ic_ext(i) >= d_ext(k+1)
              k = k + 1;
        end
        interp_mass(i) = cum_mass(k) + ...
                         atom_ext(k) * (ic_ext(i) > d_ext(k)) + ...
                         rho_ext(k) * (ic_ext(i) - d_ext(k));
    end

    % density on Eulerian mesh
    return_data = diff(interp_mass) ./ diff(ic_ext);

end