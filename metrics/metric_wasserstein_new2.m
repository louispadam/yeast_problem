function return_data = metric_wasserstein_new2(x,d1,d2)

    n = length(d2);

    w_at = 0;

    iter = 0;
    shift = 0;
    cont = true;

    inv_1 = inv_cum_dist(x,d1);

    while cont

        inv_2 = inv_cum_dist(x,circshift(d2,shift));
        inv_2_up = inv_cum_dist(x,circshift(d2,shift+1));
        %inv_2_lo = inv_cum_dist(x,circshift(d2,shift-1));
    
        w_at = lp_integrate(x,abs(inv_1-inv_2),1);
        w_up = lp_integrate(x,abs(inv_1-inv_2_up),1);
        %w_lo = lp_integrate(x,abs(inv_1-inv_2,lo),1);

        if floor(n/(2^(iter+2))) < 1
            cont = false;
        end

        if w_up - w_at <= 0
            shift = shift + floor(n/(2^(iter+2)));
        else
            shift = shift - floor(n/(2^(iter+2)));
        end

        iter = iter + 1;

    end

    return_data = [w_at, shift/n];
end

function return_data = inv_cum_dist(x,y)

    %****************************
    % Construct CDF
    %****************************
    dist_y = zeros([1,length(x)]);
    for k = 2:length(y)
        dist_y(k) = lp_integrate(x(1:k),y(1:k),1);
    end
    dist_y(1) = 0;
    dist_y(end) = 1;

    %****************************
    % Construct Inverse
    %****************************
    inv_dist_y = zeros([1,length(x)]);
    inv_x = linspace(0,1,length(x));   % same discretization as x_data
    for k = 1:length(inv_x)
        inv_dist_y(k) = x(find(dist_y>=inv_x(k),1));
    end

    return_data = inv_dist_y;

end