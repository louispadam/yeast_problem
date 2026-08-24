function W1 = metric_wasserstein1_cont(d1,d2)

    ld = length(d1);
    x = linspace(0,1-1/ld,ld);

    % Ensure column vectors
    %x  = x(:);
    %d1 = d1(:);
    %d2 = d2(:);

    N = length(x);

    % Basic checks
    %if length(d1) ~= N || length(d2) ~= N
    %    error('x, d1, and d2 must have the same length.');
    %end

    % Uniform periodic grid
    dx = x(2)-x(1);

    %if any(abs(diff(x)-dx) > 100*eps(max(1,max(abs(x)))))
    %    error('x must be uniformly spaced.');
    %end

    % Check that x represents [0,1)
    %if abs(x(1)) > 100*eps || abs(x(end) - (1-dx)) > 100*eps
    %    error('x must be a uniform grid on [0,1).');
    %end

    % --------------------------------------------------
    % Construct G = F1 - F2
    % --------------------------------------------------

    h = d1-d2;

    G = dx*cumsum([0, (h(1:end-1)+h(2:end))/2, h(end)+h(1)/2]);

    % Check equal total mass
    %total_mass_difference = G(end) + wrap;

    %if abs(total_mass_difference) > 1e-10
    %    warning('The two densities do not have equal total mass.');
    %end

    % --------------------------------------------------
    % Find median of G
    % --------------------------------------------------

    Gsort = sort(G);

    med = Gsort(ceil(N/2));

    % --------------------------------------------------
    % Compute W1
    % --------------------------------------------------

    W1 = dx*sum(abs(G-med));

end