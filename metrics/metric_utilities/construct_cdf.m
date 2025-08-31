function [return_dom, return_data] = construct_cdf(y_data,d)
%
%last updated 08/30/25 by Adam Petrucci

    % construct_cdf is not ready
    disp('WARNING: construct_cdf is not ready')

    max_val = max(y_data,'all');
    n = length(y_data);
    steps = max_val/d;
    cdf = zeros([1,steps+1]);
    for k = 0:(steps+1)
        cdf(k) = sum(y_data<k*d)/n;
    end
    return_dom = 0:d:1;
    return_data = cdf;

end