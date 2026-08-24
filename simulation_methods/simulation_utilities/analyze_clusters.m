function [return_time, return_final] = ...
    analyze_clusters(time_vec,cluster_data)

    [~,max_cluster] = max(cluster_data,[],2);
    final_cluster = max_cluster(end);
    transition_index = length(max_cluster) + 1 - ...
            find(flip(max_cluster)~=final_cluster,1,'first');
    return_time = time_vec(transition_index);
    return_final = final_cluster;

end