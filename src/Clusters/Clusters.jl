module Clusters

import LINCEGE:
    Vertices.AbstractVertices,
    Vertices.ExpansionVertices,
    Lattices.AbstractLattice,
    _NI

abstract type AbstractCluster end

expansion_vertices(cluster::AbstractCluster) = _NI("expansion_vertices")
lattice_constant(cluster::AbstractCluster) = _NI("lattice_constant")
ghash(cluster::AbstractCluster) = _NI("ghash")
cluster(old_cluster::AbstractCluster, lattice::AbstractLattice, to::Type{<:AbstractCluster}) = _NI("cluster")

function subgraphs(cluster::AbstractCluster, lattice::AbstractLattice)
    max_order = length(cluster) - 1
    roots = expansion_vertices(cluster)
    visited = Set{Subgraph}([Subgraph(ExpansionVertices(root)) for root in roots])

    function try_mark(v)
        already = v in visited
        if !already
            push!(visited, v)
        end

        if length(v) == max_order
            return false
        end
        !already
    end

    function dfs(v)
        if !try_mark(v)
            return
        end

        nbrs = neighbors(v, cluster, lattice)
        for u in nbrs
            dfs(u)
        end
    end

    for r in roots
        dfs(r)
    end

    visited

end

Base.length(cluster::AbstractCluster) = length(expansion_vertices(cluster))
Base.hash(cluster::AbstractCluster, h::Unsigned) = hash(ghash(cluster), h)
Base.isequal(g1::AbstractCluster, g2::AbstractCluster) = (ghash(g1) == ghash(g2))
Base.show(io::IO, cluster::AbstractCluster) = print(io, "Cluster with $(length(cluster)) vertices.")

include("./DirectionClusters.jl")
include("./BondClusters.jl")

export AbstractCluster,
    DirectionCluster,
    BondCluster,
    expansion_vertices,
    multiplicity,
    ghash,
    cluster
end
