struct SymmetricHasher <: AbstractHasher
        trans_hasher::TranslationHasher
        permutations::Vector{Vector{Int64}}
        connections::Union{<:AbstractConnections,Nothing}
end

function SymmetricHasher(lattice::AbstractInfiniteLattice, lattice_symmetries::Vector{Matrix{Float64}})
        trans_hasher = TranslationHasher(lattice)
        all_coords = get_coordinates(lattice)
        permutations = get_permutations(all_coords, lattice_symmetries)

        SymmetricHasher(
                trans_hasher,
                permutations,
                nothing
        )

end

function ghash(h::SymmetricHasher, lvs::LatticeVertices)
        all_hashes = Set()
        for perm in h.permutations
                push!(all_hashes, ghash(h.trans_hasher, LatticeVertices(perm[lvs])))
        end
        hash(all_hashes)
end

function ghash(h::SymmetricHasher, lvs::ExpansionVertices)
        all_hashes = Set()
        for perm in h.permutations
                push!(all_hashes, ghash(h.trans_hasher, LatticeVertices(perm[lvs])))
        end
        hash(all_hashes)
end
