function return_data = ct_sharp(d)
    return_data = @(start,stop) @(x) (start<x).*(x<stop);
end