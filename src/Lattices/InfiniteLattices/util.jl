function generate_cartesian_coordinates(dimension::Int, half_side_length::Int)
    # Forces the lattice to have a strict center point
    diameter = 2 * half_side_length + 1
    # Total number of coordinates for the entire lattice
    max_coords = diameter^dimension
    coords = repeat(0:(max_coords-1), 1, dimension)'

    for dim = 0:(dimension-1)
        coords[dim+1, :] =
            div.(coords[dim+1, :], diameter^dim) .% diameter .- half_side_length
    end

    coords
end

function generate_coordinates(max_order::Int, num_basis_elements::Int, dimension::Int)
    primitive_coordinates = generate_cartesian_coordinates(dimension, max_order)

    hcat([vcat(primitive_coordinates, repeat([i], size(primitive_coordinates, 2))') for i in 1:num_basis_elements]...)
end

function generate_weighted_adj_matrix(coordinates::AbstractMatrix{<:Int}, unit_cell::UnitCell)
    adj_matrix = zeros(Int, size(coordinates, 2), size(coordinates, 2))
    for (index, coord) in enumerate(eachcol(coordinates))
        for (neighbor, bond) in zip(find_possible_neighbors(unit_cell, coord), unit_cell.bonds)
            neighbor_index = findfirst(col -> col == neighbor, eachcol(coordinates))
            if neighbor_index !== nothing
                adj_matrix[index, neighbor_index] = bond.bond_type
                adj_matrix[neighbor_index, index] = bond.bond_type
            end
        end
    end
    adj_matrix
end

function adj_mat_to_adj_list(adj_matrix::AbstractMatrix{<:Int})
    n_coords = size(adj_matrix, 1)
    adj_list = Vector{ExpansionVertices}(undef, n_coords)
    for i in 1:n_coords
        neighbors = ExpansionVertices()
        for j in 1:n_coords
            if adj_matrix[i, j] > 0
                neighbors = union(neighbors, ExpansionVertices(j))
            end
        end
        adj_list[i] = neighbors
    end
    adj_list
end

function generate_neighbor_list(coordinates::AbstractMatrix{Int}, unit_cell::UnitCell)
    adj_matrix = generate_weighted_adj_matrix(coordinates, unit_cell)
    adj_mat_to_adj_list(adj_matrix)
end

function find_centers(coordinates::AbstractMatrix{Int})
    center_indices = findall(col -> all(x -> x == 0, col[1:end-1]), eachcol(coordinates))
    ExpansionVertices(center_indices)
end
