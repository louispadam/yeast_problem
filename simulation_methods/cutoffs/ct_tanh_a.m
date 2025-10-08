function return_data = ct_tanh_a(d)
%CT_TANH_A sets the framework for defining a bump function
%It is based on tanh, has bounded support, and is smooth.
%
%last updated 10/08/25 by Adam Petrucci
arguments (Input)
    d       % length on which to 'bump'
end

    % rescale 'bump' is on approximately length d
    dd = d/2;

    return_data = @(start,stop) @(x) ...
        0*(x < start-dd) + ...
        (start-dd < x).*(x < start+dd).*(0.5*(1+tanh((x-start)./(dd^2-(x-start).^2)))) + ...
        1*(start+dd < x).*(x< stop-dd) + ...
        (stop-dd < x).*(x < stop+dd).*(0.5*(1+tanh((stop-x)./(dd^2-(stop-x).^2)))) + ...
        0*(stop+dd < x);
end