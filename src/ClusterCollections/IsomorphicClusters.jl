struct IsomorphicClusters{H<:IsomorphicHash, C<:IsomorphicCluster} <: AbstractClusters{H, C}
    clusters::Dict{H,C}
end

IsomorphicClusters() = IsomorphicClusters(Dict{IsomorphicHash, IsomorphicCluster}())

function IsomorphicClusters(translation_clusters::TranslationClusters, lattice::AbstractLattice)
    partitioned = IsomorphicClusters()
    rlock = ReentrantLock()

    @info "Starting parallel partitioning with $(nthreads()) threads."

    @threads for (ghash, translation_cluster) in translation_clusters
        cluster = IsomorphicCluster(translation_cluster, lattice)

        lock(rlock) do
            refined[ghash(cluster)] = cluster
        end
    end

    @info "Partitioning complete with $(length(partitioned)) unique clusters."

    partitioned
end

_clusters(cs::IsomorphicClusters) = cs.clusters
Base.length(cs::IsomorphicClusters) = length(_clusters(cs))
Base.iterate(cs::IsomorphicClusters) = iterate(_clusters(cs))
Base.iterate(cs::IsomorphicClusters, state) = iterate(_clusters(cs), state)
Base.getindex(cs::IsomorphicClusters, ghash::IsomorphicHash) = getindex(cs.clusters, ghash)
Base.haskey(cs::IsomorphicClusters, ghash::IsomorphicHash) = haskey(cs.clusters, ghash)
function Base.setindex!(cs::IsomorphicClusters, ghash::IsomorphicHash, cluster::IsomorphicCluster)
    if ghash in cs
        merge!(cs.clusters[ghash], cluster)
    else
        cs.clusters[ghash] = cluster
    end
    cs
end
