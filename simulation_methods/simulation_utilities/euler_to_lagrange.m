function return_data = euler_to_lagrange(eul_coord, density)
%EULER_TO_LAGRANGE converts eulerian coordinates to lagrange coordinates.
%
%last updated 10/16/26 by Adam Petrucci
arguments (Input)
    eul_coord
    density
end
arguments (Output)
    return_data
end

    lag_mass = ([density(2:end) + density(1:end-1), density(1) + density(end)] .* ...
                [eul_coord(2:end) - eul_coord(1:end-1), mod(eul_coord(1) - eul_coord(end),1)])/2;

    return_data = [eul_coord;lag_mass];

end