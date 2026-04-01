# Utility functions that will be used in multiple modules
function pairwise_direction(coordinates::AbstractMatrix{<:Real})
        T = eltype(coordinates)
        dim, n_coords = size(coordinates)
        pw_dir = Array{T}(undef, n_coords, n_coords, dim)
        for i in 1:n_coords
                for j in 1:n_coords
                        for d in 1:dim
                                @inbounds pw_dir[i, j, d] = coordinates[d, j] - coordinates[d, i]
                        end
                end
        end
        pw_dir
end

all_lattice_symmetries = Dict([:Square => Vector{Matrix{Float64}}([[1 0; 0 1], [0 1; -1 0], [-1 0; 0 -1], [0 -1; 1 0], [0 1; 1 0], [0 -1; -1 0], [1 0; 0 -1], [-1 0; 0 1]])])
