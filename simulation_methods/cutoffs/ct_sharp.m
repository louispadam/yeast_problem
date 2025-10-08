function return_data = ct_sharp(d)
%CT_TANH_K sets the framework for defining a bump function
%This is a step function, so discontinuous
%
%last updated 10/08/25 by Adam Petrucci
arguments (Input)
    d       % not used here, only included for consistency with approximate cutoffs
end
    return_data = @(start,stop) @(x) (start<x).*(x<stop);
end