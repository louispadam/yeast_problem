function return_data = generate_parameter_space(L1, L2, L3, options)
%GENERATE_PARAMETER_SPACE helps generate a parameter space to explore with
%big_experiment. L1, L2, and L3 are vectors each containing the range of
%some parameter value to explore. This function automatically builds a
%matrix of admissible parameter sets from those vectors.
%
%last updated 09/04/26 by Adam Petrucci

arguments (Input)
    L1         % range of distances between s1 and s2
    L2         % range of distances between s1 and r1
    L3         % range of distances between r1 and r2
end
arguments (Input)
    options.delta = 0.1   % minimum gap between r2 and s1
end
arguments (Output)
    return_data    % array of parameter sets (:,4)
end

    % Extract delta
    delta = options.delta;

    % Compute tolerance value for delta to avoid floating point errors
    L_min = min([diff(L1),diff(L2),diff(L3)]/2);

    % Set up storage object
    param_space = [];

    for ind1 = 1:length(L1)

        % Compute admissible l2 given l1
        L12 = L2(L2 + min(L3) + L1(ind1) < 1 - (delta - L_min));

        for ind2 = 1:length(L12)

            % Compute admissible l3 given l1 and l2
            L123 = L3(L3 + L12(ind2) + L1(ind1) < 1 - (delta - L_min));

            for ind3 = 1:length(L123)

                % Convert gap information into parameter information
                param_space(size(param_space,1)+1,:) = ...
                    [0, L1(ind1), L1(ind1)+L2(ind2),...
                        L1(ind1)+L2(ind2)+L3(ind3)];

            end % of L3 loop

        end % of L2 loop

    end % of L1 loop

    return_data = param_space;

end