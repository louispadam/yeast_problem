function return_data = gaussian_compare(x,y,animate)
%GAUSSIAN_COMPARE Helper function for comparing the state at a given time
% to a gaussian assuming it is already a unimodal, compactly-supported
% signal 
%
%last updated 11/10/25 by Adam Petrucci

    % if desired, set up environment for animation
    if animate
        figure(4)
        clf
        hold on
    end

    % Find location of maximum value 
    [M, m] = max(y);

    % Compute center of domain
    d = round(length(y)/2);

    % Shift data to be centered on domain (otherwise periodic boundaries
    % make regression more difficult)
    if d > m
        y = [y(end-(abs(m-d)-1):end) y(1:end-abs(m-d))];
    end
    if m > d
        y = [y(abs(m-d)+1:end) y(1:abs(m-d))];
    end

    % Find region of positivity
    be = find(y>exp(-10),1,"first");
    en = find(y>exp(-10),1,"last");

    % Restrict to positive region and take logarithmic. In the case of a
    % gaussian, this should yield a quadratic
    f = y(be:en);
    lf = log(f);

    % Built-in Matlab polynomial fit
    [pcf, rrr] = polyfit(x(be:en),lf,2);

    % if desired, plot comparison of data to fit
    if animate
        plot(x(be:en),lf,LineWidth=4,DisplayName="Data")
        plot(x(be:en),polyval(pcf,x(be:en)),LineWidth=2,DisplayName="Fit")
        legend();
    end

    % return standard-deviation fitted gaussian, as well as r-squared value
    return_data = [1/sqrt(-pcf(1)), rrr.rsquared];

end
