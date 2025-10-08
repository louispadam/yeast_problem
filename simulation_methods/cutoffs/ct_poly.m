function return_data = ct_poly(d,options)
%CT_POLY sets the framework for defining a bump function
%It is polynomial, has bounded support, and is N-times continuously differentiable.
%
%last updated 10/08/25 by Adam Petrucci
arguments (Input)
    d               % length on which to 'bump'
end
arguments (Input)
    options.N = 1   % differentiability
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