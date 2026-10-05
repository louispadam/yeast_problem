function [return_x,return_cont,return_sing] = split_sing_cong(mesh,mass)
%SPLIT_SING_CONT splits a periodic Langrangian mesh with associated mass
%distribution into singular and continuous components
%
%last updated 10/05/26 by Adam Petrucci
arguments (Input)
    mesh       % Lagrangian mesh
    mass       % masses associated to left boundaries
end
arguments (Output)
    return_x      % consolidated mesh (no repeats)
    return_cont   % continuous component
    return_sing   % singular components
end

    % Set up objects
    n = length(mesh);
    return_x    = zeros(1,n);
    return_cont  = zeros(1,n);
    return_sing = zeros(1,n);

    i = 1;  % (original) mesh counter
    k = 0;  % (new) mesh counter

    while i <= n

        % increment counter
        k = k+1;

        % find string of repeats
        j = i;
        while j < n && mesh(j+1) == mesh(i)
            j = j+1;
        end

        % update mesh
        return_x(k) = mesh(i);

        % update singular component
        return_sing(k) = sum(mass(i:j-1));

        % update cont component
        if j < n
            width = mesh(j+1)-mesh(j);
        else
            width = mesh(1)+1-mesh(j);
        end
        return_cont(k) = mass(j)/width;

        % increment counter
        i = j+1;

    end

    % return output objects
    return_x    = return_x(1:k);
    return_cont = return_cont(1:k);
    return_sing = return_sing(1:k);

end