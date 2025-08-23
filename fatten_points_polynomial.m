function return_data = fatten_points_polynomial(dom,avg,std)

    dom_pts = length(dom);
    dom_len = dom(end)-dom(1);
    dom = cat(1,[ dom - dom_len , dom , dom + dom_len ]);
    total = zeros([1,3*dom_pts]);

    for mean = avg
        bump = ((dom-mean-std).*(dom-mean+std)).^2;
        cutoff = (mean-std<dom).*(dom<mean+std);
        total = total + cutoff.*bump;
    end

    total = total(1:dom_pts) + total(dom_pts+1:2*dom_pts) + total(2*dom_pts+1:end);
    return_data = total/(length(avg)*(std^5)*(16/15));

end