function return_data = ct_poly(d,options)
arguments (Input)
    d
end
arguments (Input)
    options.N = 1
end

    n = options.N;

    b = zeros([1,2*(n+1)]);
    b(1) = 1;
    f = vander_interp([d/2,-d/2],b);

    return_data = @(start,stop) @(x) ...
        0*(x < start-d/2) + ...
        (start-d/2 < x).*(x < start+d/2).*f(x-start) + ...
        1*(start+d/2 < x).*(x < stop-d/2) + ...
        (stop-d/2 < x).*(x < stop+d/2).*f(stop-x) + ...
        0*(stop+d/2 < x);

end