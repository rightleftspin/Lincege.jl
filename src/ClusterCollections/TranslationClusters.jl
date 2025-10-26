struct TranslationClusters{H<:TranslationHash, C<:TranslationCluster} <: AbstractClusters{H, C}
    clusters::Dict{H,C}
end

function TranslationClusters()
    TranslationClusters(Dict{TranslationHash,TranslationCluster}())
end

function TranslationClusters(centers::AbstractVertices, lattice::AbstractSiteExpansionLattice)
    roots = TranslationClusters{TranslationHash, TranslationCluster}()
    for v in centers
        cluster = TranslationCluster(typeof(centers)(v), lattice)
        roots[ghash(cluster)] = cluster
    end
    roots
end

"""
    TranslationClusters(lattice::AbstractSiteExpansionLattice, max_order::Int; spawn_depth::Int=2)

Parallel DFS over the subgraph tree of an arbitrary site expansion lattice.
- `lattice`: provides the neighbor clusters of each cluster
- `max_order`: maximum NLCE expansion order
- `spawn_depth`: sequentially explore until this depth; deeper levels spawn tasks
"""
function TranslationClusters(lattice::AbstractSiteExpansionLattice; spawn_depth::Int=2)
    max_order = max_order(lattice)
    roots = TranslationClusters(centers(lattice), lattice)
    visited = TranslationClusters{TranslationHash, TranslationCluster}()
    vlock = ReentrantLock()

    @info "Starting parallel DFS with $(nthreads()) threads."

    function try_mark(cluster::TranslationCluster)
        lock(vlock)
        already = cluster in visited
        if !already
            visited[ghash(cluster)] = cluster
        end
        unlock(vlock)

        if length(cluster) == max_order
            return false
        end
        !already
    end

    function dfs(cluster::TranslationCluster, depth)
        if !try_mark(cluster)
            return
        end

        nbrs = neighbor_clusters(cluster, lattice)

        if depth < spawn_depth
            for sc in nbrs
                dfs(sc, depth + 1)
            end
        else
            tasks = Task[]
            first = true
            for sc in nbrs
                if first
                    dfs(sc, depth + 1)
                    first = false
                else
                    push!(tasks, @spawn dfs(sc, depth + 1))
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

"""
    TranslationClusters(lattice::AbstractClusterExpansionLattice, max_order::Int; spawn_depth::Int=2)

Parallel DFS over the subgraph tree of an arbitrary cluster expansion lattice.
- `lattice`: provides the neighbor clusters of each cluster
- `max_order`: maximum NLCE expansion order
- `spawn_depth`: sequentially explore until this depth; deeper levels spawn tasks
"""
function TranslationClusters(lattice::AbstractClusterExpansionLattice; spawn_depth::Int=2)
    _NI("TranslationClusters for AbstractClusterExpansionLattice not yet implemented")
end

_clusters(cs::TranslationClusters) = cs.clusters
Base.length(cs::TranslationClusters) = length(_clusters(cs))
Base.iterate(cs::TranslationClusters) = iterate(_clusters(cs))
Base.iterate(cs::TranslationClusters, state) = iterate(_clusters(cs), state)
Base.getindex(cs::TranslationClusters, ghash::TranslationHash) = getindex(cs.clusters, ghash)
Base.haskey(cs::TranslationClusters, ghash::TranslationHash) = haskey(cs.clusters, ghash)
Base.setindex!(cs::TranslationClusters, ghash::TranslationHash, cluster::TranslationCluster) = setindex!(cs.clusters, ghash, cluster)
