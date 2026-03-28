struct TranslationHasher <: AbstractHasher
        hashing_matrix::Matrix{Int}
        connections::Union{<:AbstractConnections,Nothing}
end

function TranslationHasher(lattice::AbstractInfiniteLattice, connections::Union{<:AbstractConnections,Nothing})
        pwd_matrix = pairwise_direction(get_coordinates(lattice))
        hashing_matrix, max_dir = unique_direction_indices(pwd_matrix, bond_matrix(lattice))
        diag_matrix = diagm(get_labels(lattice) .+ max_dir)
        TranslationHasher(2 .^ (hashing_matrix + diag_matrix), connections)
end

TranslationHasher(lattice::AbstractInfiniteLattice) = TranslationHasher(lattice, nothing)
TranslationHasher(lattice::AbstractClusterExpansionLattice) = TranslationHasher(lattice, connections(lattice))

ghash(h::TranslationHasher, lvs::LatticeVertices) = hash(sum(h.hashing_matrix[lvs, lvs], dims=2))
ghash(h::TranslationHasher, evs::ExpansionVertices) = ghash(h, union(LatticeVertices(), h.connections[evs]))

