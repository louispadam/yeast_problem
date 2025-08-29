function return_data = lp_integrate(x_data,y_data,q)
%LP_INTEGRATE integrates using the trapezoidal method
arguments
    x_data (1,:)    % discretization of domain
    y_data (1,:)    % function values
    q               % L^q to calculate
end

    % Collect Inputs
    x = x_data;
    y = y_data;

    % Integrate!
    yp = y.*q;
    integrate = (x(2)-x(1))*(sum(yp)-0.5*yp(1)-0.5*yp(end));
    return_data = integrate^(1/q);

end