function W1 = metric_wasserstein1_lag(mesh1,mass1,mesh2,mass2)
%METRIC_WASSERSTEIN1_LAG
%
% Circular Wasserstein-1 distance between two piecewise-constant
% densities on possibly different, nonuniform periodic meshes.
%
% INPUT:
%   mesh1(i) = left boundary of cell i for measure 1
%   mass1(i) = mass in [mesh1(i),mesh1(i+1))
%
%   mesh2(i) = left boundary of cell i for measure 2
%   mass2(i) = mass in [mesh2(i),mesh2(i+1))
%
% The final cell wraps periodically through x = 1 = 0.
%
% Each density is assumed constant on each cell:
%
%       rho_i = mass_i / cell_width_i.
%
% The meshes must be strictly increasing in [0,1).
%
% Uses the circular formula
%
%       W1 = min_alpha int_0^1 |G(x)-alpha| dx,
%
% where G = F1-F2.
%
% For piecewise-constant densities, G is piecewise linear.
% The minimizing alpha is computed exactly (up to floating-point
% arithmetic) by a level sweep.


    % ============================================================
    % Put inputs in row-vector form
    % ============================================================

    mesh1 = mesh1(:).';
    mesh2 = mesh2(:).';

    mass1 = mass1(:).';
    mass2 = mass2(:).';


    % ============================================================
    % Basic checks
    % ============================================================

    n1 = length(mesh1);
    n2 = length(mesh2);

    %if length(mass1) ~= n1
    %    error('mesh1 and mass1 must have the same length.')
    %end

    %if length(mass2) ~= n2
    %    error('mesh2 and mass2 must have the same length.')
    %end

    %if any(diff(mesh1) <= 0) || any(diff(mesh2) <= 0)
    %    error('Meshes must be strictly increasing.')
    %end

    %if any(mesh1 < 0) || any(mesh1 >= 1) || ...
    %   any(mesh2 < 0) || any(mesh2 >= 1)

    %    error('Mesh points must lie in [0,1).')
    %end


    % Wasserstein distance requires equal total mass
    M1 = sum(mass1);
    M2 = sum(mass2);

    %if abs(M1-M2) > 1e-12*max([1,abs(M1),abs(M2)])
    %    error('The two measures must have equal total mass.')
    %end


    % ============================================================
    % Cell densities
    % ============================================================

    width1 = [diff(mesh1), mesh1(1)+1-mesh1(end)];
    width2 = [diff(mesh2), mesh2(1)+1-mesh2(end)];

    rho1 = mass1 ./ width1;
    rho2 = mass2 ./ width2;


    % ============================================================
    % Determine which cell of each mesh contains x = 0.
    %
    % If mesh(1)=0, cell 1 starts at zero.
    % Otherwise zero lies in the periodic final cell.
    % ============================================================

    if mesh1(1) == 0

        cell1 = 1;
        next1 = 2;

    else

        cell1 = n1;
        next1 = 1;

    end


    if mesh2(1) == 0

        cell2 = 1;
        next2 = 2;

    else

        cell2 = n2;
        next2 = 1;

    end


    % ============================================================
    % Simultaneously:
    %
    %   1. merge the two meshes,
    %   2. determine rho1-rho2,
    %   3. integrate to construct G = F1-F2.
    %
    % G(k) is the value at the left endpoint of merged segment k.
    %
    % seg_length(k) is the length of that segment.
    %
    % Because both meshes are monotone, each pointer only advances.
    % ============================================================

    max_segments = n1+n2+1;

    seg_length = zeros(1,max_segments);
    G = zeros(1,max_segments+1);

    x = 0;
    k = 0;

    while x < 1

        % Next boundary from mesh 1
        if next1 <= n1
            b1 = mesh1(next1);
        else
            b1 = 1;
        end

        % Next boundary from mesh 2
        if next2 <= n2
            b2 = mesh2(next2);
        else
            b2 = 1;
        end


        % Next point of the merged mesh
        x_new = min(b1,b2);


        % New merged segment
        k = k+1;

        seg_length(k) = x_new-x;


        % rho1-rho2 is constant on this entire segment, so
        % G is exactly linear.
        rho_diff = rho1(cell1)-rho2(cell2);

        G(k+1) = G(k) ...
               + rho_diff*seg_length(k);


        % If x_new is a mesh-1 boundary, enter the next cell
        if next1 <= n1 && x_new == b1

            cell1 = next1;
            next1 = next1+1;

        end


        % If x_new is a mesh-2 boundary, enter the next cell
        if next2 <= n2 && x_new == b2

            cell2 = next2;
            next2 = next2+1;

        end


        x = x_new;

    end


    seg_length = seg_length(1:k);
    G = G(1:k+1);


    % Equal total mass implies G(1)=G(0).
    % Keep this as a useful diagnostic.
    %if abs(G(end)) > 1e-10
    %    error('Cumulative mass difference does not close periodically.')
    %end


    % ============================================================
    % EXACT CONTINUOUS MEDIAN OF G
    % ============================================================
    %
    % We need alpha satisfying
    %
    %       |{x : G(x) <= alpha}| >= 1/2
    %
    % and
    %
    %       |{x : G(x) >= alpha}| >= 1/2.
    %
    % On a nonconstant segment with endpoint values g0,g1,
    %
    %       measure{G <= alpha}
    %
    % increases linearly from 0 to the segment length as alpha
    % moves from min(g0,g1) to max(g0,g1).
    %
    % Thus its derivative with respect to alpha is
    %
    %       segment_length / |g1-g0|
    %
    % between those two levels.
    %
    % A constant-G segment instead causes a jump in the distribution
    % of G at that level.
    %
    % We collect these level events and sweep upward in alpha.
    % ============================================================

    nseg = length(seg_length);

    % Each nonconstant segment contributes two slope-change events.
    % Each constant segment contributes one jump event.
    level        = zeros(1,2*nseg);
    slope_change = zeros(1,2*nseg);
    mass_jump    = zeros(1,2*nseg);

    nevent = 0;


    for j = 1:nseg

        g0 = G(j);
        g1 = G(j+1);

        h = seg_length(j);


        if g0 == g1

            % G is constant on an interval of length h.
            % The distribution of G therefore has an atom of
            % spatial mass h at this level.

            nevent = nevent+1;

            level(nevent) = g0;
            mass_jump(nevent) = h;


        else

            lo = min(g0,g1);
            hi = max(g0,g1);

            rate = h/(hi-lo);


            % At lo, this segment begins contributing slope
            nevent = nevent+1;

            level(nevent) = lo;
            slope_change(nevent) = rate;


            % At hi, this segment stops contributing slope
            nevent = nevent+1;

            level(nevent) = hi;
            slope_change(nevent) = -rate;

        end

    end


    level        = level(1:nevent);
    slope_change = slope_change(1:nevent);
    mass_jump    = mass_jump(1:nevent);


    % Sort level events
    [level,idx] = sort(level);

    slope_change = slope_change(idx);
    mass_jump    = mass_jump(idx);


    % ============================================================
    % Sweep through the level variable alpha.
    %
    % M(alpha) = measure{x : G(x) <= alpha}.
    %
    % Between consecutive level events M is linear.
    % We stop when M reaches 1/2.
    % ============================================================

    target = 0.5;

    M = 0;
    slope = 0;

    prev_level = level(1);

    j = 1;

    found = false;


    while j <= nevent

        this_level = level(j);


        % --------------------------------------------------------
        % Advance M continuously from prev_level to this_level
        % --------------------------------------------------------

        M_next = M ...
               + slope*(this_level-prev_level);


        if M_next >= target

            % Median lies between the two event levels.
            alpha = prev_level ...
                  + (target-M)/slope;

            found = true;
            break

        end


        M = M_next;


        % --------------------------------------------------------
        % Combine every event at this same G-level
        % --------------------------------------------------------

        total_slope_change = 0;
        total_jump = 0;

        while j <= nevent && level(j) == this_level

            total_slope_change = ...
                total_slope_change + slope_change(j);

            total_jump = ...
                total_jump + mass_jump(j);

            j = j+1;

        end


        % --------------------------------------------------------
        % Constant-G intervals create a jump in M.
        %
        % If 1/2 falls inside this jump, this level itself is
        % an exact median.
        % --------------------------------------------------------

        if M + total_jump >= target

            alpha = this_level;

            found = true;
            break

        end


        M = M + total_jump;


        % The new slope applies ABOVE this level
        slope = slope + total_slope_change;

        prev_level = this_level;

    end


    if ~found
        % Should only be needed because of floating-point arithmetic.
        alpha = max(G);
    end


    % ============================================================
    % Compute integral |G-alpha| exactly.
    %
    % G-alpha is linear on each merged segment.
    % ============================================================

    W1 = 0;


    for j = 1:nseg

        a = G(j)   - alpha;
        b = G(j+1) - alpha;

        h = seg_length(j);


        if a*b >= 0

            % No sign change.
            %
            % |G-alpha| is linear, so the trapezoidal rule
            % is exact.

            W1 = W1 ...
               + 0.5*h*(abs(a)+abs(b));


        else

            % G-alpha crosses zero inside this segment.
            %
            % The integral is the sum of two triangle areas.

            aa = abs(a);
            bb = abs(b);

            W1 = W1 ...
               + 0.5*h*(aa^2+bb^2)/(aa+bb);

        end

    end

end