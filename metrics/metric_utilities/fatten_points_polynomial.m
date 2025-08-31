function return_data = fatten_points_polynomial(domain,averages,std_dev)
%FATTEN_POINTS_POLYNOMIAL produces a continuous distribution from an
%empirical one using C^1 polynomial bumps
%
%last updated 08/30/25 by Adam Petrucci
arguments (Input)
    domain (1,:)    % domain on which to fatten
    averages (1,:)  % points to fatten
    std_dev         % thickness of bumps
end

    % Collect Inputs
    dom = domain;
    avg = averages;
    std = std_dev;

    % Construct helper objects
    dom_pts = length(dom);
    dom_len = dom(end)-dom(1);
    dom = cat(1,[ dom - dom_len , dom , dom + dom_len ]);
    total = zeros([1,3*dom_pts]);

    % Fit polynomial bump around each point
    for mean = avg
        bump = ((dom-mean-std).*(dom-mean+std)).^2;
        cutoff = (mean-std<dom).*(dom<mean+std);
        total = total + cutoff.*bump;
    end

    % Correct for periodic domain
    total = total(1:dom_pts) + total(dom_pts+1:2*dom_pts) + total(2*dom_pts+1:end);
    
    return_data = total/(length(avg)*(std^5)*(16/15));

end