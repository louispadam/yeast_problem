function [return_time, return_data, return_clock] = temp_func(initial,parameters,options)
%PARTICLE_PROOF simulates the yeast NODE using a scheme inspired from our
%proof of the mean-field limit.
%
%last updated 07/22/26 by Adam Petrucci
arguments (Input)
    initial (1,:)       % initial conditions
    parameters struct   % parameters for simulation
end
arguments (Input)
    options.Timestep = 2
end
arguments (Output)
    return_time (1,:)   % discretized time axis of simulation
    return_data (:,:)   % simulation results: [time,data]
    return_clock
end

    % Begin timer
    tic

    %****************************
    % Collect Inputs
    %****************************

    ic = initial;
    params = parameters;
    fxd_dt = options.Timestep;

    % Define Temporal parameters
    %dt = params.dt;
    t_final = params.tfin;
    msz=params.m_sz;
    ud=params.update;

    %****************************
    % Set up iteration
    %****************************

    % change coordinate system so that r2_tilde = 1 = 0
    r1_tilde = mod(params.r1 - params.r2,1);
    s1_tilde = mod(params.s1 - params.r2,1);
    s2_tilde = mod(params.s2 - params.r2,1);
    ic = mod(ic-params.r2,1);
    N = length(ic); % number of particles

    [d, labels] = sort(ic);   % iteration vector for data

    tt = 0;        % iteration value for time
    ref_lab = 0;   % label for reference particle
    ref_pos = 0;   % initial position of reference particle
    iterate_counter = 1;
    if ~any((d > s1_tilde) .* (d < s2_tilde))
        if ~any(d < s1_tilde)
            ref_lab = N;
            ref_pos = d(end);
            tt = s1_tilde - ref_pos + 1;
        else
            [val, ref_lab] = find(d < s1_tilde,'last');
            ref_pos = val;
            tt = s1_tilde - ref_pos;
        end
        d = mod(d + tt,1);
        j = find(diff(d)<0,1);
        if ~isempty(j)
            d = [d(j+1:end), d(1:j)];
            labels = [labels(j+1:end), labels(1:j)];
        end
        %iterate_counter = 2;
    else
        [ref_pos, ref_lab] = find(d > s1_tilde,1,'first');
    end

    if fxd_dt > s1_tilde
        fprintf('Using default timestep of %.4f\n',s1_tilde)
        %compute delta
        %in new coords, s1-r2=s1
        fxd_dt = s1_tilde; % otherwise known as |delta|
    end

    % Define stepping for iteration
    %steps = round(t_final/min(dt, s1_tilde) + 1);
    %sz = steps;
    %keep = 1;   % frequency with which to store iteration

    % If default time-vector is longer than permitted, replace with max
    %if sz > msz
    %    sz = msz;
    %    keep = steps/msz;
    %end

    % Define time and space vectors to store
    %time = zeros([1,sz]);
    %data = zeros([sz,length(ic)]);
    %data(1,:) = ic;

    time = zeros([1,1000]);
    time(iterate_counter) = tt;
    data = zeros([1000,length(ic)]);
    data(iterate_counter,labels) = d;

    kk = 1; % revolution counter
    %kk = 2;      % counter for storing iteration
    %here = round(keep*kk);

    % Prepare frequency of updates (if desired)
    %pb = round(linspace(2,steps,20));
    %n_pb = 1;

    %****************************
    % Iterate!
    %****************************

    % updates if desired
    %if ud
    %    fprintf("Began Simulation\n");
    %end

    new_rev = false;
    converged = false;

    %for step = 2:steps
    while ~converged && (tt < t_final)

        dt = fxd_dt;
        check_ref = ref_pos - d(mod(ref_lab - labels(end),N)+1);
        if (check_ref > 0) && check_ref < fxd_dt && ...
                              ~(abs(check_ref) < 1e-10)
           dt = check_ref;
           new_rev = true;
        end

        % determine particles that interact with S
        s_set = d(d > s1_tilde - dt & d < s2_tilde);
        % determine particles that enter S
        enter_s = s1_tilde - s_set;
        enter_s = enter_s(enter_s > 0);

        % determine particles that leave S
        exit_s  = s2_tilde - s_set;
        exit_s  = exit_s(exit_s < dt);

        % order event times and respective changes in S population
        times = [enter_s, exit_s];
        jumps = [ones(size(enter_s)),-ones(size(exit_s))];
        [times,idx] = sort(times);
        jumps = jumps(idx);

        % compute population after each event
        N0 = sum(d >= s1_tilde & d < s2_tilde);
        Ns = N0 + [0, cumsum(jumps)];

        % compute speed after each event
        speeds = 1 - params.alph*Ns/N;

        % add boundary times and compute net displacement
        tau = [0, times, dt];
        Htau = [0, cumsum(speeds .* diff(tau))];

        % determine particles that interact with R
        % in new coords, r2=1 so x<r2 is trivial
        r_array = d > r1_tilde-dt;
        r_set = d(r_array);
        r_ind = find(r_array);

        % particles' starting points relative to R
        base = max(r_set, r1_tilde);

        % compute entry times
        entry_times = zeros(size(r_set));
        enter_r = r_set < r1_tilde;
        entry_times(enter_r) = r1_tilde - r_set(enter_r);

        % compute correction for H
        Hcorr = zeros([1,length(r_set)]);
        k = 1;
        for i = flip(1:length(entry_times))
            while k < length(tau) && tau(k+1) <= entry_times(i)
                k = k + 1;
            end
            Hcorr(i) = Htau(k) + speeds(k)*(entry_times(i) - tau(k)); % H at entry time
        end

        % compute displacement of particles interacting with R
        targets = 1 - base + Hcorr;
        k = 1;
        for i = flip(1:length(r_set))
            while k < length(Htau) && Htau(k+1) <= targets(i)
                k = k + 1;
            end
            if k == length(Htau)
                d(r_ind(i)) = base(i) + Htau(end) - Hcorr(i);
            else
                exit_time = (targets(i)-Htau(k))/(speeds(k)) + tau(k);
                % technically, the next line is r2+... but r2=1=0
                d(r_ind(i)) = 1 + dt - exit_time;
            end
        end

        % compute displacement of all other particles
        d(~r_array) = d(~r_array) + dt;

        d = mod(d,1);

        j = find(diff(d)<0,1);

        if ~isempty(j)
            d = [d(j+1:end), d(1:j)];
            labels = [labels(j+1:end), labels(1:j)];
        end

        tt = tt + dt;

        % Store result at previously calculated frequency
        %if step == here
        %    time(kk) = tt;
        %    data(kk,labels) = d;
        %    kk = kk+1;
        %    here = round(keep*kk);
        %end
        if new_rev
            kk = kk+1;
            time(kk) = tt;
            data(kk,labels) = d;
            converged = metric_wasserstein1(data(kk,:),data(kk-1,:)) < 1e-4;
            disp(metric_wasserstein1(data(kk,:),data(kk-1,:)))

            %test_frame = figure(16);
            %clf(test_frame);
            %ax = axes(test_frame);
            %frame(linspace(0,1,2^11),parameters,ax,...
            %"Data",{data(kk,:)},...
            %"Meta",{struct('name',"Proof Scheme", ...
            %                 'discrete',true, ....
            %                 'color',[253,179,56]/255, ...
            %                 'thickness',parameters.del)}, ...
            %"Legend",true,...
            %"Regions",true,...
            %"Region_labels",true,...
            %"Title","Final Conditions");
            %pause(1)
        end

        % display update if desired
        %if ud && step == pb(n_pb)
        %    fprintf('Simulation Progress: %3.0f%%\n',100*step/steps)
        %    n_pb = n_pb + 1;
        %end

        iterate_counter = iterate_counter + 1;

    end

    end_time = toc;

    % Return data
    return_time = time(1:iterate_counter);
    data = data(1:iterate_counter,:);
    return_data = mod(data + params.r2,1);
    return_clock = end_time;

    % update if desired
    if ud
        fprintf('Completed Simulation in %f seconds\n',end_time)
    end

end