function return_data = generate_parameter_space_abs(L1,L2,L3,L4,options)
%GENERATE_PARAMETER_SPACE assumes L1>=0
%
%last updated 11/16/26 by Adam Petrucci
arguments (Input)
    L1
    L2
    L3
    L4
end
arguments (Input)
    options.delta = 0.1   % minimum gap between r2 and s1
    options.tol   = 1e-8
    options.min_S = 0.09
    options.min_R = 0.09
end
arguments (Output)
    return_data    % array of parameter sets (:,4)
end

    % Extract delta
    delta = options.delta;
    tol = options.tol;
    min_S = options.min_S;
    min_R = options.min_R;

    % Set up storage object
    param_space = [];

    for ind1 = 1:length(L1)

        L2_ad = L2(L2 - L1(ind1) > max(tol,min_S));

        for ind2 = 1:length(L2_ad)

            L3_ad = L3(L3 - L2_ad(ind2) > tol);

            for ind3 = 1:length(L3_ad)

                L4_ad = L4((L4 - L3_ad(ind3) > max(tol,min_R)) & ...
                           (1 - L4 > delta - tol));

                for ind4 = 1:length(L4_ad)

                    param_space(size(param_space,1)+1,:) = ...
                        [L1(ind1), L2_ad(ind2),L3_ad(ind3),L4_ad(ind4)];

                end % of L4 loop

            end % of L3 loop

        end % of L2 loop

    end % of L1 loop

    return_data = param_space;

end