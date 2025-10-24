"""
    DirectionCluster{V<:ExpansionVertices,H<:Unsigned} <: AbstractCluster

A data structure that represents a cluster on a given lattice along with its translational hash.

# Fields
- `expansion_vertices::V`: The vertices of the cluster.
- `neighbors::V`: The vertices neighboring the cluster.
- `ghash::H`: A translational hash value
"""
struct DirectionCluster{V<:ExpansionVertices,H<:Unsigned} <: AbstractCluster
    expansion_vertices::V
    neighbors::V
    ghash::H
end

neighbors(cluster::DirectionCluster) = cluster.neighbors
lattice_constant(cluster::DirectionCluster) = 1

function DirectionCluster(
    evs::ExpansionVertices,
    lattice::AbstractLattice,
)

    DirectionCluster(
        evs,
        union(ExpansionVertices(), neighbors(lattice, evs)),
        translational_hash(lattice, evs),
    )
end

function DirectionCluster(
    evs::ExpansionVertices,
    neighbors::ExpansionVertices,
    lattice::AbstractLattice,
)
    DirectionCluster(
        evs,
        neighbors,
        translational_hash(lattice, evs),
    )
end

"""
    Constructs a list of `DirectionCluster` objects from the neighbors of the given cluster.

    # Arguments
    - `cluster::DirectionCluster`: The cluster to expand around.
    - `lattice::Lattice`: The lattice to expand onto.

    # Returns
    A vector of `DirectionCluster` objects representing the neighboring clusters.
"""
function neighbor_clusters(
    cluster::DirectionCluster,
    lattice::AbstractLattice,
)
    ns = DirectionCluster[]
    for n in neighbors(cluster)
        push!(ns, neighbor_cluster(cluster, ExpansionVertices(n), lattice))
    end
    ns
end

function neighbor_cluster(
    cluster::DirectionCluster,
    n::ExpansionVertices,
    lattice::AbstractLattice,
)
    new_evs = union(expansion_vertices(cluster), n)

    DirectionCluster(
        new_evs,
        setdiff(union(neighbors(cluster), neighbors(lattice, n)), new_evs),
        lattice,
    )
end
