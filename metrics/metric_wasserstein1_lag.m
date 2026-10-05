function return_data = metric_wasserstein1_lag(mesh1,mass1,mesh2,mass2)
%METRIC_WASSERSTEIN1_LAG computes Wasserstein distance between a pair of
%periodic Lagrangian meshes. It can handle singular components
%
%last updated 10/05/26
arguments (Input)
    mesh1           % periodic Langrangian mesh
    mass1           % masses associated to to left-boundary
    mesh2           % periodic Langrangian mesh
    mass2           % masses associated to to left-boundary
end
arguments (Output)
    return_data      % Wasserstein distance
end

    % in memory of my previous version
    if any(diff(mesh1) < 0) || any(diff(mesh2) < 0)
        error('Meshes must be nondecreasing.')
    end

    %****************************
    % Construct cumulative difference
    %****************************

    % Split into singular and continuous components
    [x1,rho1,atom1] = split_sing_cont(mesh1,mass1);
    [x2,rho2,atom2] = split_sing_cont(mesh2,mass2);

    n1 = length(x1);
    n2 = length(x2);
    max_segments = n1+n2+1;
    seg_length = zeros(1,max_segments);
    G_left     = zeros(1,max_segments);
    G_right    = zeros(1,max_segments);

    G_tot = 0;

    % set up intial densities for first mesh
    if x1(1) == 0
        den1 = rho1(1);
        G_tot = G_tot + atom1(1);
        i1 = 2;
    else
        den1 = rho1(end);
        i1 = 1;
    end

    % set up intial densities for second mesh
    if x2(1) == 0
        den2 = rho2(1);
        G_tot = G_tot - atom2(1);
        i2 = 2;
    else
        den2 = rho2(end);
        i2 = 1;
    end

    x_left = 0;  % current left endpoint
    k = 0;       % (new) mesh counter
    while x_left < 1

        % determine next endpoint
        if i1 <= n1      % next value from mesh1
            b1 = x1(i1);
        else
            b1 = 1;
        end
        if i2 <= n2      % next value from mesh2
            b2 = x2(i2);
        else
            b2 = 1;
        end
        x_right = min(b1,b2);

        % increment counter
        k = k+1;

        % segment length
        seg_length(k) = x_right-x_left;

        % update continuous component of G
        G_left(k)  = G_tot;
        G_right(k) = G_left(k) + (den1-den2)*seg_length(k);
        G_tot      = G_right(k);

        % update singular component of G
        if i1 <= n1 && b1 == x_right
            G_tot = G_tot + atom1(i1);
            den1 = rho1(i1);
            i1 = i1+1;
        end
        if i2 <= n2 && b2 == x_right
            G_tot = G_tot - atom2(i2);
            den2 = rho2(i2);
            i2 = i2+1;
        end

        % update endpoint
        x_left = x_right;

    end

    % truncate data
    seg_length = seg_length(1:k);
    G_left     = G_left(1:k);
    G_right    = G_right(1:k);

    %****************************
    % Prepare cumulative G
    %****************************

    nseg = length(seg_length);       % number of segments
    level        = zeros(1,2*nseg);  % value of G
    slope_change = zeros(1,2*nseg);  % changes in slope
    mass_jump    = zeros(1,2*nseg);  % jumps in G^{-1}

    nevent = 0;   % counter for slope change
    for j = 1:nseg

        g0 = G_left(j);
        g1 = G_right(j);
        h = seg_length(j);

        if g0 == g1  % ie G constant on this segment

            nevent = nevent+1;
            level(nevent) = g0;
            mass_jump(nevent) = h;

        else

            lo = min(g0,g1);
            hi = max(g0,g1);
            rate = h/(hi-lo);

            nevent = nevent+1;           % "when G is lo, it has increasing
            level(nevent) = lo;          % clope contributed by this
            slope_change(nevent) = rate; % segment"

            nevent = nevent+1;
            level(nevent) = hi;
            slope_change(nevent) = -rate;

        end

    end

    % truncate data objects
    level        = level(1:nevent);
    slope_change = slope_change(1:nevent);
    mass_jump    = mass_jump(1:nevent);

    % Sort cumulation functions
    [level,idx]  = sort(level);
    slope_change = slope_change(idx);
    mass_jump    = mass_jump(idx);

    %****************************
    % Find median of G
    %****************************

    target = 0.5;
    M = 0;
    slope = 0;

    prev_level = level(1);

    not_found = true;
    j = 1;
    while j <= nevent && not_found

        % update current accumulation
        curr_level = level(j);
        M_next = M + slope*(curr_level-prev_level);

        % check if target achieved by slope
        if M_next >= target
            alpha = prev_level + (target-M)/slope;
            not_found = false;
        end
        M = M_next;

        % repare to update slope
        total_slope_change = 0;
        total_jump = 0;

        % sum all changes in slope at current level
        while j <= nevent && level(j) == curr_level
            total_slope_change = total_slope_change + slope_change(j);
            total_jump = total_jump + mass_jump(j);
            j = j+1;
        end

        % check if target achieved by jump
        if M + total_jump >= target && not_found
            alpha = curr_level;
            not_found = false;
        end
        M = M + total_jump;

        %nupdate slope
        slope = slope + total_slope_change;

        % update level
        prev_level = curr_level;

    end

    % ChatGPT said float-point arithmetic could cause an error here, though
    % I don't really see how
    if not_found
        alpha = level(end);
    end

    %****************************
    % Compute |G-alpha|
    %****************************

    return_data = 0;

    for j = 1:nseg

        a = G_left(j)  - alpha;
        b = G_right(j) - alpha;

        h = seg_length(j);

        if a*b >= 0 % no root implies trapezoidal sum
            return_data = return_data + 0.5*h*(abs(a)+abs(b));
        else % root implies triangular sum
            aa = abs(a);
            bb = abs(b);
            return_data = return_data + 0.5*h*(aa^2+bb^2)/(aa+bb);
        end

    end

end