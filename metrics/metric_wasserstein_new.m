function return_data = metric_wasserstein(x1,y1,x2,y2,parameters,q)
%
%last updated 11/15/25 by Adam Petrucci
arguments (Input)
    x1              % discretization of first distribution
    y1              % values of first distribution
    x2              % discretization of second distribution
    y2              % values of second distribution
    parameters
    q double        % power of metric
end

    update = parameters.update;

    % I don't think setting format is necessary
    %format long

    %load lambda lambda
    lambda = q;         % order of cost function
    %load n0 n0
    n0 = length(x1); % discretization of first distribution
    %load n1 n1
    n1 = length(x2); % discretization of second distribution

    %load X X
    %load Y Y
    X = x1;         % renaming
    Y = x2;         % renaming

    %load MX MX
    %load MY MY
    MX = y1;        % renaming
    MY = y2;        % renaming

    CX = cumsum(MX);        % cumulative distribution function of first
    CY = cumsum(MY);        % cumulative distribution function of second

    Lm = 10;        % I think this has to be calculated
    Lp = 10;        % I think this has to be calculated
    L = max(Lp,Lm);

    epsi = 1e-15;   % error tolerance

    tm = -1;        % lower theta bound
    tp =  1;        % upper theta bound
    tc = (tm+tp)/2; % midpoint theta

    pasfini = 1;   % boolean for while loop

    % checking helper functions work
    test = 0;
    if test
        Ctest1 = C(0  ,X,Y,CX,CY,n0,n1,lambda);
        Ctest2 = C(0.2,X,Y,CX,CY,n0,n1,lambda);
        s = rand();
        tau	= 1e-4;
        [Dd,Dd]	= dC(s,X,Y,CX,CY,n0,n1,lambda);
        Ddapprox = (C(s+tau,X,Y,CX,CY,n0,n1,lambda) ...
                      - C(s,X,Y,CX,CY,n0,n1,lambda)) / tau;
        Dgapprox = (C(s,X,Y,CX,CY,n0,n1,lambda) ...
                  - C(s-tau,X,Y,CX,CY,n0,n1,lambda)) / tau;
    end

    iter = 0;

    % if desired, display update
    if update
        fprintf(2,'iteration 0 | tm = %f | tp = %f | tc = %f | cout = %f\n',...
                tm,tp,tc,C(tc,X,Y,CX,CY,n0,n1,lambda));
    end

    while pasfini

        iter = iter + 1;

        [dCp,dCm] = dC(tc,X,Y,CX,CY,n0,n1,lambda); % calculate derivatives
        pasfini	= (dCp*dCm)>0; % check for (no) extrema

        % check if not extrema but tolerance reached
        if ((tp-tm)<epsi/L) && pasfini

            % if desired, warn tolerance reached
            if update
                fprintf('\n ***** tp-tm petit ***** \n ')
            end

            % calculate derivatives at theta bounds
            [dCptp,dCmtp] = dC(tp,X,Y,CX,CY,n0,n1,lambda);
            [dCptm,dCmtm] = dC(tm,X,Y,CX,CY,n0,n1,lambda);

            % caculate cost at theta bounds
            Ctp	= C(tp,X,Y,CX,CY,n0,n1,lambda);
            Ctm	= C(tm,X,Y,CX,CY,n0,n1,lambda);

            % if not already within tolerance, recalculate theta
            if abs(dCptm-dCmtp)>epsi % I changed 0.001 to epsi
                tc = (Ctp-Ctm + dCptm*tm - dCmtp*tp)/(dCptm-dCmtp);
            end

            pasfini = 0;    % end loop

        % check if not extrema but tolerance not reached
        elseif pasfini > 0

            % if 0+ derivative negative, shift lower theta bound up
            if dCp<0
                tm = tc;
            else % if 0+ derivative positive, shift upper theta bound down
                tp = tc;
            end

            % update theta
            tc = .5*(tm+tp);

        % if already found extrema then nothing to do
        % I don't think we need this last line
        else
            tc = tc;
        end

        % if desired, display update
        if update
            fprintf(2,'iteration %i | tm = %f | tp = %f | tc = %f | dCp = %f | dCm =%f | cout = %f\n', ...
                    iter,tm,tp,tc,dCp,dCm,C(tc,X,Y,CX,CY,n0,n1,lambda))
        end
        
        % I think these were just print statements for debugging
        %cout=C(tc,X,Y,CX,CY,n0,n1,lambda)
        %tc

    end
end


function valC = C(t,X,Y,CX,CY,n0,n1,lambda)
% Calculate left and right derivatives of cost function

    Ip = CY-(t-floor(t))>=0;
    In = CY-(t-floor(t))< 0;

    valF0 = ones(n0+n1,n0)*diag(CX);
    valF1t = ones(n0+n1,n1)*diag([CY(Ip)-(t-floor(t)),CY(In)-(t-floor(t))+1]);
    vsort = [0,sort([valF0(1,:),valF1t(1,:)])];

    Yt = [ Y(Ip)+floor(t) ,  Y(In)+1+floor(t) ];
    Yt = [ Yt, Yt(1)+1];

    vk0 = diag(.5*(vsort(2:n0+n1+1)+vsort(1:n0+n1)))*ones(n0+n1,n0);
    vk1 = diag(.5*(vsort(2:n0+n1+1)+vsort(1:n0+n1)))*ones(n0+n1,n1);

    [vxk,nxk] = max(vk0< valF0 ,[] , 2  );
    [vyk,nyk] = max([vk1,vk1(:,1)]<[valF1t,ones(n0+n1,1)] ,[] , 2);
    xk = reshape(X (nxk), n0+n1 , 1 );
    yk = reshape(Yt(nyk), n0+n1 , 1 );
    
    % I don't know what these were for
    %size(xk);
    %size(yk);
    
    valC = [vsort(2:n0+n1+1)-vsort(1:n0+n1)]*cout(xk,yk,lambda);
    
end

function [valdcp,valdcm] = dC(t,X,Y,CX,CY,n0,n1,lambda)
% Calculate cost

    Ip = CY-(t-floor(t))>=0;
    In = CY-(t-floor(t))< 0;

    valF0 = ones(n1,n0)*diag(CX);
    valF1t = diag([CY(Ip)-(t-floor(t)),CY(In)-(t-floor(t))+1])*ones(n1,n0);
    vsort = [0,sort([valF0(1,:),[valF1t(:,1)]'])];
    diff = .5*(abs(vsort(2:n0+n1+1) - vsort(1:n0+n1)));
    epsilon	= min(diff(diff>0));

    X = [ X , X(1)+1];
    Yt = [ Y(Ip)+floor(t) , Y(In)+1+floor(t) ];
    Yt = [ Yt, Yt(1)+1];

    [vxk,nxk] = max(valF1t <= valF0, [] , 2);
    [vxkm,nxkm] = max([valF1t,valF1t(:,1)] <= ([valF0-epsilon,ones(n1,1)]) , [] , 2);
    xk = X (nxk' );
    xkm	= X (nxkm');

    valdcp = sum(cout(xk ,Yt(2:n1+1),lambda) - cout(xk ,Yt(1:n1),lambda));
    valdcm	= sum(cout(xkm,Yt(2:n1+1),lambda) - cout(xkm,Yt(1:n1),lambda));

end

function res=cout(x,y,lambda)
% the cost metric

    res = (abs(x-y)).^lambda;

end