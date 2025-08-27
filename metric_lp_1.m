function return_data = metric_lp_1(dis_y,cont_x,cont_y,q)

    poly_thic = linspace(0,0.3,201);
    poly_thic = poly_thic(2:end);

    trans = linspace(0,1,201);
    trans = trans(1:end-1);

    solns = zeros([length(poly_thic),length(trans)]);

    for k1 = 1:length(poly_thic)
        pt = poly_thic(k1);
        dis_thic = fatten_points_polynomial(cont_x,dis_y,pt);
        for k2 = 1:length(trans)
            t = trans(k2);
            act_trans = round(t*length(cont_x));
            cont_trans = cat(2,cont_y(act_trans+1:end),cont_y(1:act_trans));
            diff = abs(dis_thic-cont_trans);
            norm = lp_integrate(cont_x,diff,q);
            solns(k1,k2) = norm;
        end
    end

    inf_norm = min(solns,[],'all');
    [t1,t2] = find(solns == inf_norm);
    best_thic = poly_thic(t1);
    best_trans = trans(t2);

    return_data = [inf_norm,best_thic,best_trans];

end