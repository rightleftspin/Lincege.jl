module ClusterSets

import LINCEGE:
    _NI,
    Vertices.AbstractVertices,
    Lattices.AbstractLattice,
    Clusters.AbstractCluster,
    Clusters.DirectionCluster

abstract type AbstractClusterSet{C<:AbstractCluster} end

hashes(cs::AbstractClusterSet) = _NI("hashes")
get_order(cs::AbstractClusterSet, order::Int) = _NI("get_order")
Base.push!(cs::AbstractClusterSet{C}, cluster::C) where C <:AbstractCluster = _NI("Base.push!")

Base.length(cs::AbstractClusterSet) = length(hashes(clusters))

Base.in(cs::AbstractClusterSet, cluster::AbstractCluster) = ghash(cluster) in cs
Base.in(cs::AbstractClusterSet, ghash::Unsigned) = ghash in hashes(cs)

Base.haskey(cs::AbstractClusterSet, cluster::AbstractCluster) = haskey(cs, ghash(cluster))
Base.haskey(cs::AbstractClusterSet, ghash::Unsigned) = haskey(hashes(cs), ghash)

Base.show(io::IO, cs::AbstractClusterSet) = print(io, "ClusterSet with $(length(cs)) clusters.")

export AbstractClusterSet,
    hashes,
    get_order

include("BondClusterSets.jl")
include("DirectionClusterSets.jl")
end
