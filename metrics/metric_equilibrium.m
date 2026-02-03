function return_data = metric_equilibrium(x,t_data,p_data,c_data)

    d_p_d = histcounts(p_data(end,:), x);
    d_p_d = cat(2,d_p_d,[0]);
    lp = lp_integrate(x,d_p_d,1);
    d_p_d = d_p_d/lp;

    steps_in_rev = sum(t_data>t_data(end)-1);

    trans_dists = zeros([1,steps_in_rev]);

    for i = 1:steps_in_rev
         trans_dists(i) = lp_integrate(x,abs(c_data(end-i+1,:)-d_p_d),1);
    end

    [min_val, min_ind] = min(trans_dists);

    return_data = [min_val,min_ind];

end