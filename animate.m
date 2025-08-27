function return_data = animate(x,y,params,names,pn)

    dt=params.dt; % time step
    tmax = params.tfin;

    names(1)
    names(2)

    tt=0;
    ptfac=params.fr;   % 1/animation rate

    pp = params.pr;

    while (tt<tmax)

        if (abs((fix(ptfac*tt)-ptfac*tt))/ptfac<dt)

            ind = floor(tt/dt)+1;

            frame(x,y,params,names,pn,ind,tt);

            pause(pp)

        end

        tt=tt+dt;

    end

    return_data = 1;

end