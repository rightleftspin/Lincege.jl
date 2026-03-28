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

function unique_direction_indices(pw_dir::AbstractArray{<:Real,3}, bond_matrix::AbstractMatrix{Int})
        n_coords = size(pw_dir, 1)
        dir_to_index = Dict{Vector{Float64},Int}()
        index_matrix = zeros(Int, n_coords, n_coords)

        for i in 1:n_coords, j in 1:n_coords
                if i != j && bond_matrix[i, j] != 0
                        dir = round.(pw_dir[i, j, :], digits=8)
                        if !haskey(dir_to_index, dir)
                                dir_to_index[dir] = length(dir_to_index) + 1
                        end
                        index_matrix[i, j] = dir_to_index[dir]
                end
        end

        return index_matrix, length(dir_to_index)
end

function get_permutations(coords::Matrix{Float64}, syms::Vector{Matrix{Float64}})
        n_sites = size(coords, 2)

        permutations = Vector{Vector{Int64}}()
        for sym in syms
                perm = zeros(Int64, n_sites)
                transformed = sym * coords
                for i in 1:n_sites
                        j = findfirst(k -> isapprox(transformed[:, i], coords[:, k]), 1:n_sites)
                        isnothing(j) || (perm[i] = j)
                end
                push!(permutations, perm)
        end

        permutations
end

