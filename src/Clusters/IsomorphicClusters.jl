struct IsomorphicCluster{V<:AbstractVertices, H<:IsomorphicHash, P<:AbstractPermutation} <: AbstractCluster{V, H}
    vertices::V
    init_perm::P
    ghash::H
    permutations::Vector{P}
end

function IsomorphicCluster(translation_cluster::TranslationCluster, lattice::AbstractLattice)
    vs = vertices(translation_cluster)
    ghash = IsomorphicHash(lattice, vs)
    IsomorphicCluster(vs, EmptyPermutation(), ghash, EmptyPermutation[])
end

vertices(cluster::IsomorphicCluster) = cluster.vertices
ghash(cluster::IsomorphicCluster) = cluster.ghash
permutations(cluster::IsomorphicCluster) = cluster.permutations
lattice_constant(cluster::IsomorphicCluster) = length(permutations(cluster)) + 1
merge!(cluster1::IsomorphicCluster, cluster2::IsomorphicCluster) = append!(permutations(cluster1), permutations(cluster2))
