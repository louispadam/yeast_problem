function return_data = lp_integrate(x,y,q)

    yp = y.*q;
    integrate = (x(2)-x(1))*(sum(yp)-0.5*yp(1)-0.5*yp(end));
    return_data = integrate^(1/q);

end