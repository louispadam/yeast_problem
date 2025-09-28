function return_data = sample(x_data,dist,num_points)
%Sampling Sample a given probability distribution to construct an
%approximating empirical measure.
%
%This differs from construct_em in that I implement built-in Matlab
%functions, guided by ChatGPT
%
%last updated 09/24/25 by Adam Petrucci
arguments (Input)
    x_data (1,:) double % discretization of domain for y_data
    dist                % distribution (function) to be sampled from
    num_points          % number of points in sample
end

    % Collect required inputs (for notational convenience)
    n = num_points;

    % Construct CDF of dist
    %cumtrapz is cumulative integration via trapezoid method
    pdf_vals = dist(x_data);
    cdf_vals = cumtrapz(x_data, pdf_vals);

    % Invert CDF and sample
    %interp1 is linear interpolation; note switching (x,y) in the above
    %simulates the inverse
    p_unif = rand([1,n]);
    p_dist = interp1(cdf_vals, x_data, p_unif);

    % Return data
    return_data = p_dist;
end