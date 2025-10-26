module ClusterSets

using Base.Threads

import LINCEGE:
    _NI,
    GraphHashes.AbstractGraphHash,
    GraphHashes.IsomorphicHash,
    GraphHashes.TranslationHash,
    GraphHashes.VertexHash,
    Vertices.AbstractVertices,
    Lattices.AbstractLattice,
    Lattices.AbstractSiteExpansionLattice,
    Lattices.AbstractClusterExpansionLattice,
    Clusters.AbstractCluster,
    Clusters.TranslationCluster,
    Clusters.IsomorphicCluster,
    Clusters.Subgraph

abstract type AbstractClusters{H<:AbstractGraphHash, C<:AbstractCluster} end

Base.length(cs::AbstractClusters) = _NI("Base.length")
Base.iterate(cs::AbstractClusters) = _NI("Base.iterate")
Base.iterate(cs::AbstractClusters, state) = _NI("Base.iterate")
Base.getindex(cs::AbstractClusters, ghash::AbstractGraphHash) = _NI("Base.getindex")
Base.setindex!(cs::AbstractClusters, ghash::AbstractGraphHash, cluster::AbstractCluster) = _NI("Base.setindex!")
Base.haskey(cs::AbstractClusters, ghash::AbstractGraphHash) = _NI("Base.haskey")

Base.in(cs::AbstractClusters, ghash::AbstractGraphHash) = haskey(cs, ghash)
Base.in(cs::AbstractClusters, cluster::AbstractCluster) = haskey(cs, ghash(cluster))
Base.show(io::IO, cs::AbstractClusters) = print(io, "Cluster Collection with $(length(cs)) clusters.")

function get_order(cs::AbstractClusters, order::Int)
    result = typeof(cs)()
    for (ghash, cluster) in cs
        if length(cluster) == order
            result[ghash] = cluster
        end
    end
    result
end

include("TranslationClusters.jl")
include("IsomorphicClusters.jl")
include("Subgraphs.jl")

export AbstractClusters, TranslationClusters, IsomorphicClusters, Subgraphs, get_order

end
