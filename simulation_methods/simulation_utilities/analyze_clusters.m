function [return_time, return_final] = ...
    analyze_clusters(time_vec,cluster_data)
%ANALYZE_CLUSTERS helps analyze cluster data. Given the an array of cluster
%modes, it determines the dominating mode (the final number of clusters)
%and the time at which this mode becomes the dominant.
%
%last updated 08/25/26 by Adam Petrucci
arguments (Input)
    time_vec       % time data
    cluster_data   % magnitude of cluster modes at each time step
end
arguments (Output)
    return_time    % time of transition to final number of clusters
    return_final   % number of clusters at end of experiment
end

    [~,max_cluster] = max(cluster_data,[],2);
    final_cluster = max_cluster(end);
    transition_index = length(max_cluster) + 1 - ...
            find(flip(max_cluster)~=final_cluster,1,'first');
    return_time = time_vec(transition_index);
    return_final = final_cluster;

end