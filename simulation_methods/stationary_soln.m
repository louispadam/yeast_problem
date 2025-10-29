function return_data = stationary_soln(x,params)

    r1 = params.r1;
    r2 = params.r2;
    s1 = params.s1;
    s2 = params.s2;
    alph = params.alph;
    eps = params.eps;
    
    r = r2-r1;
    S = s2-s1;

    % Calculate inital guess in eps\to\infty case
    a0 = (sqrt((1+alph*S)^2+4*alph*S*(r-1))-(1+alph*S))/(2*alph*S*(r-1));
    g0 = 1-alph*a0*S;

    s1i = mod(s1 - r1,1);
    s2i = mod(s2 - r1,1);

    % Run built-in Matlab solver
    for_solve = @(v) fun(v,eps,alph,r,s1i,s2i);

    %options = optimoptions('fsolve','Display','none','PlotFcn',@optimplotfirstorderopt);
    options = optimoptions('fsolve','Display','none');
    soln = fsolve(for_solve,[a0,g0],options);

    au = soln(1);
    gu = soln(2);
    bu = ss_b(au,gu,r);
    cu = ss_c(au,gu,r);

    % Construct solution
    ff = @(x) (au + bu*exp((r-x)/eps^2)).*(x>=r) + ...
              (au/gu + cu*exp(-gu*x/eps^2)).*(x<r);

    yy = ff(x);
    rs = length(x)-round(r1*length(x));
    yy1 = yy(1:rs);
    yy2 = yy(rs+1:end);
    yy = [yy2 yy1];

    return_data = yy;

end

% Calculate b from a, g, and r
function return_data = ss_b(a,g,r)
    num = a*(1/g-1)*(1-exp(-g*r/eps^2));
    den = (1-exp((r*(1-g)-1)/eps^2));
    return_data = num/den;
end

% Calculate c from a, g, and r
function return_data = ss_c(a,g,r)
    num = a*(1/g-1)*(exp((r-1)/eps^2)-1);
    den = 1-exp((r*(1-g)-1)/eps^2);
    return_data = num/den;
end

% Calculate a from b, c, g, and r
function return_data = ss_a(b,c,g,r)
    num = c-b*exp((r-1)/eps^2);
    den = 1-1/g;
    return_data = num/den;
end

% Calculate g from a, b, and parameters
function return_data = ss_g(a,b,alpha,s1,s2)
    S = s2-s1;
    return_data = 1-alpha*(a*S-b*eps^2*(exp(-s2/eps^2)-eps(-s1/eps^2)));
end

% Map to find zero of
function f = fun(v,eps,alpha,r,s1,s2)

    ai = v(1);
    gi = v(2);

    bi = ss_b(ai,gi,r);
    ci = ss_c(ai,gi,r);

    f(1) = ai*(1-r)+ai*r/gi-bi*eps^2*(exp((r-1)/eps^2)-1)-ci*eps^2*(exp(-gi*r/eps^2)-1)/gi-1;
    f(2) = ss_g(ai,bi,alpha,s1,s2) - gi;

end