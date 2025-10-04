function return_data = ct_tanh_k(d)
    return_data = @(start,stop) @(x) 0.25*(tanh(2*(x-start)/d)+1).*(tanh(2*(stop-x)/d)+1);
end