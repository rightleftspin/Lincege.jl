struct TranslationHasher <: AbstractHasher
    hashing_matrix::AbstractMatrix{Int}
end

function TranslationHasher(lattice::SiteExpansionLattice)
    pwd_matrix = pairwise_direction(get_coordinates(lattice))
    hashing_matrix::Matrix{Int}, max_dir = unique_direction_indices(pwd_matrix, bond_matrix(lattice))
    diag_matrix = diagm(get_labels(lattice) .+ max_dir)
    TranslationHasher(2 .^ (hashing_matrix + diag_matrix))
end

ghash(h::TranslationHasher, evs::ExpansionVertices) = hash(sum(h.hashing_matrix[evs, evs], dims=2))

struct IsomorphicHasher <: AbstractHasher
    hashing_matrix::AbstractMatrix{Int}
    labels::AbstractVector{Int}
end

function IsomorphicHasher(lattice::SiteExpansionLattice)
    IsomorphicHasher(bond_matrix(lattice), )
end

ghash(h::TranslationHasher, evs::ExpansionVertices) = hash(sum(h.hashing_matrix[evs, evs], dims=2))
