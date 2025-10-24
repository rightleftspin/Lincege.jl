struct BondCluster{V<:AbstractVertices} <: AbstractCluster
    expansion_vertices::V
    clusters_composed::Vector{<:AbstractCluster}
    permutation::Vector{Vector{Int}}
    ghash::Unsigned
end
