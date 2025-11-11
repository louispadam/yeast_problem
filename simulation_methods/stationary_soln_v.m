function return_data = stationary_soln_v(x,params)

    r1 = params.r1;
    r2 = params.r2;
    s1 = params.s1;
    s2 = params.s2;
    alph = params.alph;
    
    r = r2-r1;
    S = s2-s1;

    % Calculate inital guess in eps\to\infty case
    a = (sqrt((1+alph*S)^2+4*alph*S*(r-1))-(1+alph*S))/(2*alph*S*(r-1));
    g = 1-alph*S*a;

    % Construct solution
    ff = @(x) a.*(x>=r) + ...
              (a/g).*(x<r);

    yy = ff(x);
    rs = length(x)-round(r1*length(x));
    yy1 = yy(1:rs);
    yy2 = yy(rs+1:end);
    yy = [yy2 yy1];

    return_data = yy;

end