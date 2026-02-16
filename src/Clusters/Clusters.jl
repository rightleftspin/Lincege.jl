module Clusters

using Base.Threads

import LINCEGE:
    _NI,
    Vertices.ExpansionVertices,
    Lattices.centers,
    Lattices.max_order,
    Lattices.neighbors,
    Lattices.AbstractLattice,
    Lattices.SiteExpansionLattice

abstract type AbstractCluster end
abstract type AbstractHasher end
abstract type AbstractClusterSet{C<:AbstractCluster,H<:AbstractHasher} end

# Cluster Methods
Base.length(c::AbstractCluster) = _NI("Base.length")
Base.hash(c::AbstractCluster, h::UInt) = _NI("Base.hash")

Base.isequal(c1::C, c2::C) where {C<:AbstractCluster} = c1 == c2
Base.:(==)(c1::C, c2::C) where {C<:AbstractCluster} = (hash(c1) == hash(c2))

# Cluster Set Methods
Base.length(cs::AbstractClusterSet)::Int = _NI("Base.length")
Base.in(cluster::C, cs::AbstractClusterSet{C,H}) where {C<:AbstractCluster,H} = _NI("Base.in")
Base.iterate(cs::AbstractClusterSet) = _NI("Base.iterate")
Base.iterate(cs::AbstractClusterSet, state) = _NI("Base.iterate")
Base.push!(cs::AbstractClusterSet{C,H}, c::C) where {C<:AbstractCluster,H} = _NI("Base.push!")
add_cluster!(cs::AbstractClusterSet, cluster::AbstractCluster) = _NI("add_cluster!")

include("Hashers.jl")
include("Cluster.jl")

function clusters_from_lattice!(clusters::AbstractClusterSet{C,H}, lattice::SiteExpansionLattice; spawn_depth::Int=3) where {C<:AbstractCluster,H}
    max_depth = max_order(lattice)
    roots = [C(center, clusters) for center in centers(lattice)]
    vlock = ReentrantLock()

    function try_mark(cluster::AbstractCluster)
        lock(vlock)
        already = cluster in clusters
        if !already
            push!(clusters, cluster)
        end
        unlock(vlock)

        if length(cluster) == max_depth
            return false
        end
        !already
    end

    function dfs(cluster::AbstractCluster, depth)
        if !try_mark(cluster)
            return
        end

        if depth < spawn_depth
            for ev in neighbors(lattice, cluster.evs)
                dfs(C(union(cluster.evs, ev), clusters), depth + 1)
            end
        else
            tasks = Task[]
            first = true
            for ev in neighbors(lattice, cluster.evs)
                neighbor_cluster = C(union(cluster.evs, ev), clusters)
                if first
                    dfs(neighbor_cluster, depth + 1)
                    first = false
                else
                    push!(tasks, @spawn dfs(neighbor_cluster, depth + 1))
                end
            end
            for t in tasks
                fetch(t)
            end
        end
    end

    for c in roots
        dfs(c, 0)
    end

    clusters
end

function clusters_from_clusters!(new_clusters::AbstractClusterSet{C,H}, old_clusters::AbstractClusterSet) where {C<:AbstractCluster,H}
    for old_cluster in old_clusters
        add_cluster!(new_clusters, old_cluster)
    end
    new_clusters
end

end
