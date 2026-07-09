function return_data = stationary_soln_vm(x,params)
%STATIONARY_SOLN_VM generates the stationary solution in the case of a
%linear influence function with diffusion as grained by the vector x by
%solving fixed point problem (really it finds the zero)
%
%last updated 07/06/26 by Adam Petrucci
arguments (Input)
    x           % input vector on which to compute solution
    params      % system parameters
end

    r1 = params.r1;
    r2 = params.r2;
    s1 = params.s1;
    s2 = params.s2;
    alph = params.alph;
    eps = params.eps;

    % Correct parameters to compute in r1=0 setup
    r   = mod(r2 - r1,1);
    s1i = mod(s1 - r1,1);
    s2i = mod(s2 - r1,1);
    S = s2i-s1i;

    % Calculate inital guess in eps\to\infty case
    au = (sqrt((1+alph*S)^2+4*alph*S*(r-1))-(1+alph*S))/(2*alph*S*(r-1));
    gu = 1-alph*au*S;

    % The stationary solution without diffusion-correction
    ff = @(x) (au).*(x>=r) + ...
              (au/gu).*(x<r);

    if eps ~= 0

        % Run built-in Matlab solver
        for_solve = @(v) fun(v,eps,alph,r,s1i,s2i);

        %options = optimoptions('fsolve','Display','none','PlotFcn',@optimplotfirstorderopt);
        options = optimoptions('fsolve','Display','none');
        soln = fsolve(for_solve,[au,gu],options);

        au = soln(1);
        gu = soln(2);

        bu = ss_b(au,gu,r);
        cu = ss_c(au,gu,r);

        % Construct solution
        ff = @(x) ff(x) + (bu*exp((x-r)/eps^2)).*(x>=r) + ...
                          (cu*exp(gu*x/eps^2)).*(x<r);

    end

    yy = ff(x);
    rs = length(x)-round(r1*length(x));
    yy1 = yy(1:rs);
    yy2 = yy(rs+1:end);
    yy = [yy2 yy1];

    return_data = yy;

end

% Calculate b from a, g, and r
function return_data = ss_b(a,g,r)
    num = a*(1/g-1)*(1-exp(g*r/eps^2))
    den = (1-exp((r*(g-1)+1)/eps^2))
    return_data = num/den;
end

% Calculate c from a, g, and r
function return_data = ss_c(a,g,r)
    num = a*(1/g-1)*(exp((r-1)/eps^2)-1)
    den = 1-exp((r*(g-1)+1)/eps^2)
    return_data = num/den;
end

% Calculate g from a, b, and parameters
function return_data = ss_g(a,b,alpha,s1,s2)
    S = s2-s1;
    return_data = 1-alpha*(a*S-b*eps^2*(exp(-s2/eps^2)-eps(-s1/eps^2)));
end

% Map to find zero of
function f = fun(v,eps,alpha,r,s1,s2)

    v

    ai = v(1);
    gi = v(2);

    bi = ss_b(ai,gi,r)
    ci = ss_c(ai,gi,r)

    f(1) = ai*(1-r)+ai*r/gi-bi*eps^2*(exp((1-r)/eps^2)-1)-ci*eps^2*(exp(gi*r/eps^2)-1)/gi-1;
    f(2) = ss_g(ai,bi,alpha,s1,s2) - gi;

    f

end