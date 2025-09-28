function return_data = find_stationary(r,s1,s2,alpha,eps)

% Compute initial guess
a0_num = sqrt((1+alpha*(s2-s1))^2+4*alpha*(s2-s1)*(r-1))-(1+alpha*(s2-s1));
a0_den = 2*alpha*(s2-s1)*(r-1);
a0 = a0_num/a0_den;
gamma0 = 1-alpha*(s2-s1)*a0;
b0 = a0*(1/gamma0-1);
c0 = -b0;

old = [a0, b0, c0, gamma0];
new = [0, 0, 0, 0];

f = @(x) exp(x/(eps*eps));

for k = 1:10

    a = old(1);
    b = old(2);
    g = old(4);

    a_num = b*(f(1-r-g)-1);
    a_den = (1/g-1)*(f(1-r-g)+f(r-1)-2);
    anew = a_num/a_den;

    b_num = (a*(1-r)+a*r/g-1)*(f(r-1)-f(-g))-...
            eps*eps*a*(1/g-1)*(1-f(r-1))*(f(-g*r)-1)/g;
    b_den = eps*eps*(f(1-r)-1)*(f(1-r)-f(-g*r));
    bnew = b_num/b_den;

    gnew = 1-alpha*a*(s2-s1)+alpha*b*eps*eps*(f(-s2)-f(s1));

    c_num = anew*(1/gnew-1)*(1-f(r-1));
    c_den = f(1-r)-f(-gnew*r);
    cnew = c_num/c_den;

    new = [anew, bnew, cnew, gnew];

    norm(new-old)

    old = new;

end

return_data = new;

end