function [d,mass,nmerge] = check_cell_widths(d,mass,tol)
%CHECK_CELL_WIDTHS merges cells when they become too small for computer
%precision, as approximated by tol
%
%last updated 10/02/26 by Adam Petrucci
arguments (Input)
    d          % cell boundaries
    mass       % cell masses m(i) corresponds to left boundary d(i)
    tol        % width to merge
end
arguments (Output)
    d          % updated boundaries
    mass       % updated masses
    nmerge     % total number of merges
end

    % compute widths
    width = [diff(d), d(1)+1-d(end)];

    N = length(d);
    ind = 1;
    nmerge = 0;

    % loop across whole vector
    while ind < N+1

        % get current width
        curr = width(ind);

        % check width admissability
        if curr < tol

            % determine neighboring cells
            ind_left = mod(ind-2,N)+1;
            ind_right = mod(ind,N)+1;

            % determine which cell to merge with
            if width(ind_left) < width(ind_right) % merge left

                % merge current into left, delete current
                width(ind_left) = width(ind_left) + curr;
                mass(ind_left) = mass(ind_left) + mass(ind);
                mass(ind) = [];
                d(ind) = [];

            else % merge right

                % merge right into current, delete current
                width(ind) = curr + width(ind_right);
                mass(ind) = mass(ind) + mass(ind_right);
                mass(ind_right) = [];
                d(ind_right) = [];

            end

            % update vector length
            N = N-1;
            nmerge = nmerge + 1;

        else

            % if no merge is necessary, move on to next cell
            ind = ind + 1;

        end

    end

end