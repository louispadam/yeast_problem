function return_data = postprocess_big_experiment_old(exp_data)
%POSTPROCESS_BIG_EXPERIMENT_OLD
%
%last updated 10/20/26 by Adam Petrucci

    % Might be explicit formula for this, I should look into it later
    total_sets = 0;
    for ind1 = 1:length(exp_data)
        for ind2 = 1:length(exp_data{ind1})
            for ind3 = 1:length(exp_data{ind1}{ind2})
                total_sets = total_sets + 1;
            end
        end
    end

    num_nums = length(exp_data{1}{1}{1}{3});
    trials = length(exp_data{1}{1}{1}{3}{1});

    data_f = cell([1,3]);
    data_f{2} = cell([1,2]);
    data_f{3} = cell([1,2]);
    data_f{3}{1} = NaN([4,trials * total_sets * num_nums]);
    data_f{3}{2} = cell([1,num_nums]);

    for k = 1:num_nums
        data_f{3}{2}{k} = NaN([trials,length(exp_data{1}{1}{1}{3}{k}{1}{3})]);
    end

    param_array = NaN([4,total_sets]);
    mf_stats = NaN([4,total_sets]);
    mf_states = NaN([length(x),total_sets]);

    temp_check = 1;

    cc = 1;
    for ind1 = 1:length(exp_data)
        for ind2 = 1:length(exp_data{ind1})
            for ind3 = 1:length(exp_data{ind1}{ind2})

                if isempty(exp_data{ind1}{ind2}{ind3}{2}{4})
                    mf_stat = [exp_data{ind1}{ind2}{ind3}{2}{1},...
                               exp_data{ind1}{ind2}{ind3}{2}{2},...
                               -1,...
                               exp_data{ind1}{ind2}{ind3}{2}{5}];
                    temp_check = temp_check + 1;
                    if temp_check == 3
                        ind1
                        ind2
                        ind3
                    end
                else
                    mf_stat = [exp_data{ind1}{ind2}{ind3}{2}{1},...
                               exp_data{ind1}{ind2}{ind3}{2}{2},...
                               exp_data{ind1}{ind2}{ind3}{2}{4},...
                               exp_data{ind1}{ind2}{ind3}{2}{5}];
                end

                mf_state = exp_data{ind1}{ind2}{ind3}{2}{3};

                param_set = exp_data{ind1}{ind2}{ind3}{1};
                param_array(:,cc) = [param_set.s1,...
                                     param_set.s2,...
                                     param_set.r1,...
                                     param_set.r2];
                mf_stats(:,cc) = mf_stat;
                mf_states(:,cc) = mf_state;

                for n = 1:num_nums
                    for t = 1:trials

                        data_f{3}{2}{n}(t,:) = exp_data{1}{1}{1}{3}{n}{t}{3};

                        if isempty(exp_data{ind1}{ind2}{ind3}{3}{n}{t}{4})
                            node_stat = [exp_data{ind1}{ind2}{ind3}{3}{n}{t}{1},...
                                     exp_data{ind1}{ind2}{ind3}{3}{n}{t}{2},...
                                     -1,...
                                     exp_data{ind1}{ind2}{ind3}{3}{n}{t}{5}];
                        else
                            node_stat = [exp_data{ind1}{ind2}{ind3}{3}{n}{t}{1},...
                                         exp_data{ind1}{ind2}{ind3}{3}{n}{t}{2},...
                                         exp_data{ind1}{ind2}{ind3}{3}{n}{t}{4},...
                                         exp_data{ind1}{ind2}{ind3}{3}{n}{t}{5}];
                        end

                        data_f{3}{1}(:,(n-1)*total_sets*trials + ...
                                       (cc-1)*trials + t) = node_stat;

                    end
                end

                cc = cc + 1;
            end
        end
    end

    data_f{1} = param_array;
    data_f{2}{1} = mf_stats;
    data_f{2}{2} = mf_states;

    return_data = data_f;

end