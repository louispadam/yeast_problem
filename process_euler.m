function return_data = process_euler(dom,time,data,std)
    
    fe_f = zeros([1,length(time),length(dom)]);
    for k = 1:length(time)
        fe_f(1,k,:) = fatten_points_polynomial(dom,data(:,k).',std);
    end

    return_data = fe_f;
    
end