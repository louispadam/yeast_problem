function sampler = construct_sampler(density,options)
%MAKE_CIRCLE_SAMPLER Construct a Monte Carlo sampler for a distribution
%characterized by the given density function
%
% last updated 07/31/26 by Adam Petrucci
arguments (Input)
    density                     % density function for distribution
end
arguments (Input)
    options.Num_grid = 2^15     % grid points for inverse CDF
end

    x = linspace(0,1,options.Num_grid);

    pdf = density(x);

    F = cumtrapz(x,pdf);

    % Remove repeated CDF values arising from zero-density regions
    keep = [true diff(F) > 1e-12];

    x = x(keep);
    F = F(keep);

    % Exact endpoints
    x(1) = 0;
    x(end) = 1;
    F(1) = 0;
    F(end) = 1;

    % Monotone inverse-CDF interpolant
    Finv = griddedInterpolant(F,x,'pchip');

    sampler = @(n) mod(Finv(rand(1,n)),1);

end