function return_data = metric_wasserstein1(x,y)

    % don't need if assume sorted
    x = sort(mod(x(:),1));
    y = sort(mod(y(:),1));

    n = length(x);
    m = length(y);

    events = [x; y];
    jumps  = [ones(n,1)/n; -ones(m,1)/m];

    [events,~,ic] = unique(events);
    jumps = accumarray(ic,jumps);

    [events,idx] = sort(events);
    jumps = jumps(idx);

    events2 = [events; events(1)+1];

    lengths = diff(events2);

    gvals = cumsum(jumps);

    % should check that this is actually median
    [gSort,idx] = sort(gvals);
    wSort = lengths(idx);
    cw = cumsum(wSort);
    med = gSort(find(cw>=0.5,1));
    return_data = sum(lengths .* abs(gvals - med));

end