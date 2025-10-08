% function return_data = ct_exp(d)
%     f = @(t) exp(-1./((d/2)^2-t.^2));
%     n_f = integral(f,-(d/2),(d/2));
%     return_data = @(start,stop) @(x) ...
%         0*(x <= start-d/2) + ...
%         (start-d/2 < x).*(x < start+d/2).*arrayfun(@(y) integral(@(z) f(z-start), start-d/2, y) / n_f, x) + ...
%         1*(start+d/2 <= x).*(x<= stop-d/2) + ...
%         (stop-d/2 < x).*(x < stop+d/2).*arrayfun(@(y) 1-integral(@(z) f(stop-z), stop-d/2, y) / n_f, x) + ...
%         0*(stop+d/2 <= x);
% end

% this function has issues, but I think they're about machine accuracy, not
% a math error in the design. I made need to do something mildly clever to have it
% calculate

function return_data = ct_exp(d)
%CT_EXP sets the framework for defining a bump function
%It is exponential, has bounded support, and is smooth.
%It currently fails for moderate values of d due to (I suspect) issues with machine tolerance.
%
%last updated 10/08/25 by Adam Petrucci
    return_data = @(start,stop) @(x) arrayfun(@(y) helper(start,stop,d,y),x);
end

function return_data = helper(start,stop,d,x)

    dd = d/2;

    f = @(t) exp(-1./((dd)^2-t.^2));
    n_f = integral(f,-(dd),(dd));

    g = @(t) integral(@(z) f(z),-dd,t)/n_f;

    n = (x>start-dd) + (x>=start+dd) + (x>stop-dd) + (x>=stop+dd);

    % the cases need to be split up because some blow-up outside their
    % intended regimes
    switch n
        case 1
            return_data = g(x-start);
        case 2
            return_data = 1;
        case 3
            return_data = 1-g(x-stop);
        otherwise
            return_data = 0;
    end

end