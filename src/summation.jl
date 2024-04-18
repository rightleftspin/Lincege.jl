function summation_dict(cluster::Cluster, pivot_tag::Int, pivot_tag_inner::Int)
    subcluster_info = Dict()
    for subcluster in subclusters(cluster)
        subcluster_info[get_tag(subcluster, pivot_tag)] = multiplicity(subcluster)
    end

    for subcluster in subclusters(cluster)
        for sub_subcluster in subclusters(subcluster)
            subcluster_info[get_tag(sub_subcluster, pivot_tag_inner)] -= subcluster_info[get_tag(subcluster, pivot_tag)] * multiplicity(sub_subcluster)
        end
    end

    return subcluster_info
end
