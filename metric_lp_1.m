function return_data = metric_lp_1(discrete_y,continuous_x,continuous_y,q)
%METRIC_LP_1 quantifies the difference between one empirical and one
%continuous distribution via the L^p norm by thickening the empirical
%distriubtion. It takes the infimum over a range of thickness and
%translations.
arguments
    discrete_y (1,:)    % collection of particles in empirical measure
    continuous_x (1,:)  % discretization of continuous distribution
    continuous_y (1,:)  % values of continuous distribution
    q                   % power to be taken in L^q norm
end

    % Collect Inputs
    dis_y = discrete_y;
    cont_x = continuous_x;
    cont_y = continuous_y;

    % Range of thicknesses to check
    poly_thic = linspace(0,0.3,201);
    poly_thic = poly_thic(2:end);

    % Range of translations to check
    trans = linspace(0,1,201);
    trans = trans(1:end-1);

    % Object to store results
    solns = zeros([length(poly_thic),length(trans)]);

    % Calculate norms
    for k1 = 1:length(poly_thic) % choose thickness
        pt = poly_thic(k1);
        dis_thic = fatten_points_polynomial(cont_x,dis_y,pt);
        for k2 = 1:length(trans) % choose translation
            t = trans(k2);
            act_trans = round(t*length(cont_x));
            cont_trans = cat(2,cont_y(act_trans+1:end),cont_y(1:act_trans));
            diff = abs(dis_thic-cont_trans);
            norm = lp_integrate(cont_x,diff,q);
            solns(k1,k2) = norm;
        end
    end

    % Determine infimum
    inf_norm = min(solns,[],'all');
    [t1,t2] = find(solns == inf_norm);
    best_thic = poly_thic(t1);  % optimal thickness
    best_trans = trans(t2);     % optimal translation

    return_data = [inf_norm,best_thic,best_trans];

end