function [return_time, return_data, return_clock] = cont_proof(initial, params)
%CONT_PROOF simulates the yeast MF using a scheme inspired from our
%proof of the mean-field limit. It is built on the method of
%characteristics. An effort was made to mirror structure of the analogous
%NODE scheme.
%
%last updated 08/18/26 by Adam Petrucci

    % Begin timer
    tic

    % Collect inputs
    N = length(initial);
    dx = 1/N;
    x = linspace(0,1-dx,N);
    ic = initial;

    % Change coordinate system so that s1_tilde = 1 = 0.
    r1_tilde = mod(params.r1 - params.s1,1);
    r2_tilde = mod(params.r2 - params.s1,1);
    s2_tilde = mod(params.s2 - params.s1,1);

    % Shift initial conditions according to change in coordinates
    % This may incorporate some error; I should consider adjusting delta
    [~, shift_ind] = min(abs(x-params.s1));
    d = circshift(ic,-(shift_ind-1));

    tt = 0;                   % current time

    % Automatically use optimal timestep
    % Recall 1 = s1_tilde
    fxd_dt = 1-r2_tilde;

    % temp var for num of iterations while debugging
    sss = 2;

    % Set up storage objects
    time = zeros([1,sss]);
    data = zeros([sss,length(ic)]);
    data(1,:) = ic;

    % Find S in new coordinates
    s_set = find(x < s2_tilde);            % [0,s2)  indices

    r_all = find(x > r1_tilde);            % (r1,1)  indices
    r_after = find(x > r2_tilde);          % (r2,1)  indices
    r_set = setdiff(r_all,r_after);        % (r1,r2] indices
    %the_rest = setdiff(1:length(x),r_all); % [0,r1]  indices

    %[~,ind] = max(s_set);        % first index in S
    ind = length(s_set);
    %ind = ind - 1;

    x_steps = floor(fxd_dt/dx);  % spatial steps / time step (unit speed)

    % Correction matrix for quadrature over S
    % Once S(t0) is known, each mini-step effectively translates the S
    % region, so we remove conc at ends and add conc at beginning
    % to maintian trapezoidal integration
    % Most of this nonsense is accounting for issues with mod
    corr = [-1,0,ind-1,ind] - (0:(x_steps-1)).';
    %temp = corr(:,3:4) < 1;
    %corr(:,3:4) = corr(:,3:4) - temp;
    %corr(1,2) = 1;
    corr = mod(corr,length(x)) + 1;

    % Iterate!
    for m = 2:sss

        %disp(m)

        % timestep for this iteration (adjust later for convergence
        % criterion)
        dt = x_steps * dx;

        % compute signal during timestep
        N0 = trapz(d(s_set))*dx;
        Ns = N0 + [0,cumsum(0.5*sum([1,1,-1,-1].*d(corr),2)).']*dx;

        % compute speeds and displacements
        speeds = 1 - params.alph*Ns(1:end-1);
        times = 0:dx:dt;  % S-event times
        H = [0, cumsum(speeds .* diff(times))];

        % these objects pertain to R-events
        entry_times = zeros(size(r_all));
        exit_times = dt*ones(size(r_all));
        z0 = zeros(size(r_all));   % 'starting' locations

        % particles after R exit predictably (speed 1)
        % other particles can use exit_time = dt
        exit_times(r_after - min(r_all) + 1) = dt - (x(r_after) - r2_tilde);
        
        % determine displacement-function at time of exit to correct actual
        % displacement later
        Hcorr = zeros(size(r_all));
        k = length(times);
        for i = 1:length(r_all)
            while k > 1 && times(k-1) >= exit_times(i)
                k = k - 1;
            end
            Hcorr(i) = H(k-1) + speeds(k-1)*(exit_times(i) - times(k-1)); % H at exit time
        end

        % correct exit times for later. those that don't exit (ie are in R)
        % get NaN
        exit_times(r_set - min(r_all) + 1) = NaN;
        exit_times = exit_times(~isnan(exit_times));

        % determine where to compute displacement from (ic or r1) and set
        % target displacement accordingly with correction for exit time
        base = min(x(r_all), r2_tilde);
        targets = Hcorr - (base - r1_tilde);

        k = length(times);
        for i = 1:length(r_all)
            while k > 1 && H(k-1) >= targets(i)
                k = k - 1;
            end
            if k == 1 % particle didn't enter R in this step
                z0(i) = base(i) - Hcorr(i);
                entry_times(i) = NaN;
            else              % particle entered R in this step
                %entry_times(i) = times(k);
                entry_times(i) = times(k-1) + (targets(i)-H(k-1))/speeds(k-1);
                z0(i) = r1_tilde - entry_times(i);
            end
        end

        % save those that entered
        that_enter = r_all(~isnan(entry_times));  % indices

        % correct entry times for later
        entry_times = entry_times(~isnan(entry_times));

        % r_after and that_exit should be the same the same
        %d_new = circshift(d,1-shift_ind);
        %d_new = d;

        d_new = circshift(d,x_steps);

        %[~,inds_r] = min(abs(z0 - x.'));
        %d_new(r_all) = d(inds_r);
        x_ext = [x-1, x, x+1];
        d_ext = [d,   d, d];
        z0_mod = mod(z0,1);
        d_new(r_all) = interp1(x_ext,d_ext,z0_mod,'linear');

        %[~,inds_exit] = min(abs(exit_times - times.'));
        Ns_exit = interp1(times,Ns,exit_times,'linear');
        d_new(r_after)    = d_new(r_after) .* ...
                             exp(-params.alph * Ns_exit);

        %[~,inds_enter] = min(abs(entry_times - times.'));
        Ns_enter = interp1(times,Ns,entry_times,'linear');
        d_new(that_enter) = d_new(that_enter) .* ...
                              exp(params.alph * Ns_enter);
        %max(x(that_enter))

        figure(5)
        plot(x(that_enter),exp(params.alph * Ns_enter))

        %[~,inds_rest] = min(abs(mod(x(the_rest) - dt,1) - x.'));
        %d_new(the_rest) = d(inds_rest);
        
        %check = that_enter;
        %figure(8)
        %plot(check/N,ones(size(check)),...
        %      Marker='o')
        %xlim([0,1])

        %figure(9)
        %plot(times,H)
        %min(H);

        %figure(10)
        %plot(x(r_all),x(r_all)-z0,'.-')
        %hold on
        %yline(r1_tilde)
        %yline(r2_tilde)
        %xline(r1_tilde)
        %xline(r2_tilde)
        %hold off

        target_R = r2_tilde-r1_tilde;
        t_star = interp1(H,times,target_R);
        nabla = dt-t_star;
        %fprintf('t_star = %.6f\n',t_star)
        %fprintf('nabla  = %.6f\n',nabla)

        d = d_new;

        tt = tt + dt;

        time(m) = tt;
        data(m,:) = circshift(d,shift_ind - 1);

    end

    end_time = toc;

    % Return data
    return_clock = end_time;
    return_time = time;
    return_data = data;

end