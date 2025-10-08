function return_data = ct_tanh_k(d)
%CT_TANH_K sets the framework for defining a bump function
%This is the bump suggested by Keith Promislow.
%It is based on tanh, has unbounded support, and is smooth.
arguments (Input)
    d       % length on which to 'bump'
end

    % rescale 'bump' is on approximately length d
    dd = d/4;

    return_data = @(start,stop) @(x) 0.25*(tanh((x-start)/dd)+1).*(tanh((stop-x)/dd)+1);
end