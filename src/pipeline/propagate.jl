"""
This is step three of the pipeline. In this step, the algorithm takes in an array of clusters and a pruning function. Using this information, the algorithm will grow subclusters of the clusters by treating the clusters as the 'underlying cluster' in the grow step. It will then take these subclusters and prune them according to the pruning function. Lastly, it will collect the clusters into a hashmap containing the input clusters and a list of all their 
corresponding subclusters.
"""

"""
Main function in this connection. Takes clusters and finds all possible 
subclusters of each cluster, pruning them along the way.

Inputs: 
      clusters: hashmap of cluster hashes and their related
      clusters and multiplicities

      pruning_fxn: function that will be applied to the array of 
      subclusters, needs to return a hash of the subcluster and 
      a potential permutation

Output:
      Hashmap relating the cluster to a hashmap of its subclusters
      and their corresponding multiplicities
"""
function propogate(
    pruning_fxn::Function,
    clusters::AbstractDict{Integer,Tuple{<:AbstractCluster,<:Real,<:AbstractDict}},
)

    for (hash, (cluster, mult)) in clusters
        # Grow all possible subclusters smaller than the cluster then 
        # prune all subclusters and add them to the output dictionary
        clusters[hash] = (cluster, mult, prune_cluster(pruning_fxn, grow(cluster)))
    end

    clusters
end

"""
TODO: Parallel version of this algorithm
"""
function propogate_par() end

"""
Collects clusters of the same type and counts their corresponding 
multiplicities. Need a separate method since we don't want to
keep all the subclusters in memory

Inputs: 
      clusters: array of clusters to be pruned

      pruning_fxn: function that will be applied to each cluster,
      needs to return a hash of the cluster and the permutation 
      required to convert the cluster to it's isomorphic form

Output:
      Hashmap relating the hash of the cluster to the corresponding 
      multiplicity of the cluster.
"""
function prune_cluster(pruning::Function, clusters::AbstractVector{<:AbstractCluster})

    # Initialize the empty output dictionary
    cluster_mult = Dict{Integer,Real}()

    for cluster in clusters
        hash, _ = pruning(cluster)
        # Add this information to the output dictionary
        cluster_mult[hash] = get(cluster_mult, hash, 0) + 1
    end

    cluster_mult
end

