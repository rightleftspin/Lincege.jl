#struct Clusters{C<:AbstractClusters} <: AbstractClusters
#    cluster_hashes::OrderedSet{Unsigned}
#    clusters::Vector{C}
#end
#
#Base.eltype(cs::Clusters{C}) where {C<:AbstractClusters} = C
#
#hashes(cs::Clusters) = bcs.cluster_hashes
#clusters(cs::Clusters) = bcs.clusters
#
#multiplicities(cs::Clusters) = multiplicity.(clusters(cs))
#permutations(bcs::BondClusters) = permutations.(clusters(bcs))
#
#Base.in(bcs::BondClusters, cluster::BondCluster) = haskey(hashes(bcs), ghash(cluster))
#Base.in(bcs::BondClusters, ghash::Unsigned) = haskey(hashes(bcs), ghash)
#
#Base.haskey(bcs::BondClusters, cluster::BondCluster) = haskey(hashes(bcs), ghash(cluster))
#Base.haskey(bcs::BondClusters, ghash::Unsigned) = haskey(hashes(bcs), ghash)
#
#function BondClusters()
#    BondClusters(
#        OrderedSet{Unsigned}(),
#        ExpansionVertices[],
#        Real[],
#        Vector{Vector{Int}}[],
#    )
#end
#
#function add_cluster!(bcs::BondClusters, cluster::BondCluster)
#    error("add_cluster! not implemented for BondClusters")
#end
#
#function partition(old_clusters::AbstractClusterSet, lattice::AbstractLattice, final::Type{<:AbstractCluster})
#    partitioned = SimplePartition{final, eltype(old_clusters)}()
#    rlock = ReentrantLock()
#
#    @info "Starting parallel partitioning with $(nthreads()) threads."
#
#    @threads for old_cluster in old_clusters
#        cluster = cluster(old_cluster, lattice, final)
#
#        lock(rlock) do
#            push!(partitioned, cluster)
#        end
#    end
#
#    @info "Partitioning complete with $(length(partitioned)) unique clusters."
#
#    partitioned
#end
#
