"""
    Subgraph{V<:AbstractVertices} <: AbstractCluster

A data structure that represents a subgraph of a given cluster.

# Fields
- `vertices::V`: The vertices of the subgraph.
- `parent_vertices::V`: The vertices of the parent cluster.
- `neighbors::V`: The vertices neighboring the cluster.
"""
struct Subgraph{V<:AbstractVertices, H<:VertexHash} <: AbstractCluster{V, H}
    vertices::V
    ghash::H
    parent_vertices::V
    neighbors::V
end

function Subgraph(
    vs::AbstractVertices,
    cluster::AbstractCluster,
    lattice::AbstractLattice,
)

    Subgraph(
        vs,
        VertexHash(vs),
        vertices(cluster),
        intersect(neighbors(lattice, vs), vertices(cluster)),
    )
end

"""
    Constructs a list of `Subgraph` objects from the neighbors of the given cluster.

    # Arguments
    - `cluster::Subgraph`: The cluster to expand around.
    - `lattice::AbstractLattice`: The lattice to expand onto.

    # Returns
    A vector of `Subgraph` objects representing the neighboring clusters.
"""
function neighbor_subgraphs(
    subgraph::Subgraph,
    lattice::AbstractLattice,
)
    ns = Subgraph[]
    for n in neighbors(cluster)
        push!(ns, neighbor_subgraph(subgraph, ExpansionVertices(n), lattice))
    end
    ns
end

function neighbor_subgraph(
    subgraph::Subgraph,
    n::AbstractVertices,
    lattice::AbstractLattice,
)
    new_vs = union(vertices(subgraph), n)

    Subgraph(
        new_vs,
        VertexHash(new_vs),
        parent_vertices(subgraph),
        setdiff(union(neighbors(subgraph), intersect(neighbors(lattice, n), parent_vertices(subgraph))), new_vs),
    )
end

vertices(cluster::Subgraph) = cluster.vertices
ghash(cluster::Subgraph) = cluster.ghash
parent_vertices(cluster::Subgraph) = cluster.parent_vertices
neighbors(cluster::Subgraph) = cluster.neighbors
lattice_constant(cluster::Subgraph) = 1
