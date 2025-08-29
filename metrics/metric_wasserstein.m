function return_data = metric_wasserstein(dis_y,cont_x,cont_y)
    
    % metric_wasserstein is not ready
    disp('WARNING: metric_wasserstein is not ready')

    poly_thic = linspace(0,0.5,101);
    poly_thic = poly_thic(2:end);

    trans = linspace(0,1,101);
    trans = trans(1:end-1);

    inf_norm = length(dis(y))*lp_integrate(cont_x,cont_y,q);
    best_thic = -1;
    best_trans = -1;

    for pt = poly_thic
        for t = trans
            dis_thic = fatten_points_polynomial(cont_x,dis_y,pt);
            act_trans = round(t*length(cont_x)/100);
            cont_trans = cat(2,cont_y(act_trans+1:end),cont_y(1:act_trans));

            [cdf_x, cont_cdf] = construct_cdf(cont_trans,0.001);
            [cdf_x, dis_cdf] = construct_cdf(dis_thic,0.001);

            diff = abs(dis_cdf-cont_cdf);
            norm = lp_integrate(cdf_x,diff,1);

            if norm < inf_norm
                inf_norm = norm;
                best_thic = pt;
                best_trans = t;
            end
        end
    end

    return_data = [inf_norm,best_thic,best_trans];

end