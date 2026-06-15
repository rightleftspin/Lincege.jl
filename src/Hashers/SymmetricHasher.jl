struct SymmetricHasher{C<:Union{<:AbstractConnections,Nothing}} <: AbstractHasher
        trans_hasher::TranslationHasher
        permutations::Vector{Vector{Int64}}
        connections::C
end

function SymmetricHasher(lattice::AbstractInfiniteLattice, lattice_symmetries::Vector{Matrix{Float64}}, connections::Union{<:AbstractConnections,Nothing})
        trans_hasher = TranslationHasher(lattice)
        all_coords = get_coordinates(lattice)
        # `get_permutations` may produce incomplete permutations (entries set to 0)
        # when a symmetry maps a boundary site outside the finite lattice cube.
        # For infinite-lattice use this is safe: clusters are grown from the center
        # and never reach the corners, so no cluster site will land on a zero entry.
        # For finite lattices the incomplete permutations must be filtered out before
        # constructing this hasher, otherwise `perm[lvs]` below will index with 0.
        permutations = get_permutations(all_coords, lattice_symmetries)

        SymmetricHasher(
                trans_hasher,
                permutations,
                connections
        )

end

SymmetricHasher(lattice::AbstractInfiniteLattice, lattice_symmetries::Vector{Matrix{Float64}}) = SymmetricHasher(lattice, lattice_symmetries, nothing)
SymmetricHasher(lattice::AbstractClusterExpansionLattice, lattice_symmetries::Vector{Matrix{Float64}}) = SymmetricHasher(lattice, lattice_symmetries, connections(lattice))

n_unique_sites(h::SymmetricHasher) = n_unique_sites(h.trans_hasher)
function ghash(h::SymmetricHasher, lvs::LatticeVertices)
        all_hashes = Set()
        for perm in h.permutations
                push!(all_hashes, ghash(h.trans_hasher, LatticeVertices(perm[lvs])))
        end
        hash(all_hashes)
end

ghash(h::SymmetricHasher{StrongClusterConnections}, evs::ExpansionVertices) = ghash(h, h.connections[evs])

function ghash(h::SymmetricHasher{WeakClusterConnections}, evs::ExpansionVertices)
        all_hashes = Set()
        for perm in h.permutations
                lvs, mask = h.connections[evs]
                new_lvs = perm[lvs]
                hm = h.trans_hasher.hashing_matrix[sort(new_lvs), sort(new_lvs)]
                hm[mask[sortperm(new_lvs), sortperm(new_lvs)]] .= 0
                push!(all_hashes, hash(sum(hm, dims=2)))
        end

        hash(all_hashes)
end
