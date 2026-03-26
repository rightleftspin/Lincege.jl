struct SymmetricHasher <: AbstractHasher
        trans_hasher::TranslationHasher
        permutations::Vector{Vector{Int64}}
        connections::Union{<:AbstractConnections,Nothing}
end

function SymmetricHasher(lattice::SiteExpansionLattice, lattice_symmetries::Vector{Matrix{Float64}})
        trans_hasher = TranslationHasher(lattice)
        all_coords = get_coordinates(lattice)
        permutations = get_permutations(all_coords, lattice_symmetries)

        SymmetricHasher(
                trans_hasher,
                permutations,
                nothing
        )

end


function ghash(h::SymmetricHasher, evs::ExpansionVertices)
        h = if isnothing(h.connections)
                all_hashes = Set()
                for perm in h.permutations
                        push!(all_hashes, ghash(h.trans_hasher, ExpansionVertices(perm[evs])))
                end
                h = hash(all_hashes)
        else
                _NI("symmetric ghash for cluster expansions")
        end

        h
end
