function return_data = frame(x,y,params,names,pn,ind,tt)

    if (nargin<6)
        ind = 1;
        tt = 0;
    end

    %L=2*pi;
    %a0=params.s1*2*L-L;
    %a1=params.s2*2*L-L;
    %b0=params.r1*2*L-L;
    %b1=params.r2*2*L-L;
    %del=params.del;

    a0 = params.s1;
    a1 = params.s2;
    b0 = params.r1;
    b1 = params.r2;
    del = params.del;

    ChiR = 0.25*(tanh((x-b0)/del)+1).*(tanh((b1-x)/del)+1);
    ChiS = 0.25*(tanh((x-a0)/del)+1).*(tanh((a1-x)/del)+1);

    figure(pn)
    clf
    hold on

    colors = ['b','g','r','m','y'];

    si = size(y);
    for k = 1:si(1)
        u = squeeze(y(k,ind,:));
        plot(x,u,colors(k),'lineWidth',3,'DisplayName',names(k))
    end
    m = max(y(:,ind,:),[],"all");

    %si = size(y);
    %if si(2) > 1
    %    for k = 1:si(2)
    %        u = y(ind,k,:);
    %        plot(x,u,'b','linewidth',3)
    %    end
    %    m = max(1,y(ind,:,:));
    %else
    %    plot(x,y(ind,1,:),'b','linewidth',3)
    %    m = max(y);
    %end

    plot(x,ChiR,'r.','linewidth',3)
    plot(x,ChiS,'k-','linewidth',3)

    %axis([-L L 0 m.*1.1])
    axis([0 1 0 m.*1.1]);
    tlt=sprintf('Time =%2.5g',tt);
    title(tlt,'Fontsize',18)
    legend()

    return_data = 1;

end