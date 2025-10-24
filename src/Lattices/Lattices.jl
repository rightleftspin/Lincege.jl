module Lattices

using Base.Threads
import LINCEGE:
    _NI,
    Vertices.AbstractVertices

abstract type AbstractLattice end
abstract type AbstractSiteExpansionLattice <: AbstractLattice end
abstract type AbstractClusterExpansionLattice <: AbstractLattice end

centers(lattice::AbstractLattice) = _NI("centers") # Note that centers should be a subtype of AbstractClusterSet
max_order(lattice::AbstractLattice) = _NI("max_order")
min_order(lattice::AbstractLattice) = _NI("min_order")
neighbors(lattice::AbstractLattice, vs::AbstractVertices) = _NI("neighbors")

"""
    parallel_dfs(lattice::AbstractLattice, max_order::Int; spawn_depth::Int=2)

Parallel DFS over the subgraph tree of an arbitrary lattice.
- `lattice`: provides the neighbor clusters of each cluster
- `max_order`: maximum NLCE expansion order
- `spawn_depth`: sequentially explore until this depth; deeper levels spawn tasks

Returns `visited`, where:
- `visited` is an `AbstractClusterSet` (thread-safe via a lock)
"""
function parallel_dfs(lattice::AbstractLattice; spawn_depth::Int=2)
    max_order = max_order(lattice)
    roots = centers(lattice)
    visited = typeof(roots)()
    vlock = ReentrantLock()

    @info "Starting parallel DFS with $(nthreads()) threads."

    function try_mark(v)
        lock(vlock)
        already = v in visited
        if !already
            push!(visited, v)
        end
        unlock(vlock)

        if length(v) == max_order
            return false
        end
        !already
    end

    function dfs(v, depth)
        if !try_mark(v)
            return
        end

        nbrs = neighbors(lattice, v)

        if depth < spawn_depth
            for u in nbrs
                dfs(u, depth + 1)
            end
        else
            tasks = Task[]
            first = true
            for u in nbrs
                if first
                    dfs(u, depth + 1)
                    first = false
                else
                    push!(tasks, @spawn dfs(u, depth + 1))
                end
            end
            for t in tasks
                fetch(t)
            end
        end
    end

    for r in roots
        dfs(r, 0)
    end

    @info "DFS complete with $(length(visited)) clusters."

    visited
end

export AbstractLattice,
    centers,
    max_order,
    neighbors,
    parallel_dfs

include("./SiteExpansionLattices.jl")
end
