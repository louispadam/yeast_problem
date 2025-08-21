function return_data = metric_lp_1(dis_y,cont_x,cont_y,q)

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
            diff = abs(dis_thic-cont_trans);
            norm = lp_integrate(cont_x,diff,q);
            if norm < inf_norm
                inf_norm = norm;
                best_thic = pt;
                best_trans = t;
            end
        end
    end

    return_data = [inf_norm,best_thic,best_trans];
end