function return_data = process_pseudospectral(dom,time,data)
    ps_f = zeros([1,length(time),length(dom)]);
    for k = 1:length(time)
        ps_f(1,k,:) = data(k,:);
    end
    return_data = ps_f;
end