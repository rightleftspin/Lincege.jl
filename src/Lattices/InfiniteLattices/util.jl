function generate_cartesian_coordinates(dimension::Int, half_side_length::Int)
        # Forces the lattice to have a strict center point
        r = -half_side_length:half_side_length
        hcat(vec(collect.(Iterators.product(fill(r, dimension)...)))...)
end

function generate_coordinates(max_order::Int, num_basis_elements::Int, dimension::Int)
        primitive_coordinates = generate_cartesian_coordinates(dimension, max_order)

        hcat([vcat(primitive_coordinates, repeat([i], size(primitive_coordinates, 2))') for i in 1:num_basis_elements]...)
end

function generate_coord_index(coordinates::Matrix{Int})
        coord_index = Dict{Vector{Int},Int}()
        for (i, col) in enumerate(eachcol(coordinates))
                coord_index[col] = i
        end
        coord_index
end

function generate_adj_matrix(coordinates::AbstractMatrix{Int}, unit_cell::UnitCell)
        coord_index = generate_coord_index(coordinates)
        adj_matrix = zeros(Int, size(coordinates, 2), size(coordinates, 2))
        for (index, col) in enumerate(eachcol(coordinates))
                for bond in unit_cell.bonds
                        if bond.site1 == col[end]
                                neighbor_coord = neighbor_site(bond, col)
                                ni = get(coord_index, neighbor_coord, nothing)
                                if ni !== nothing
                                        adj_matrix[index, ni] = bond.bond_type
                                        adj_matrix[ni, index] = bond.bond_type
                                end
                        end
                end
        end

        adj_matrix
end

function generate_adj_matrix(coordinates::AbstractMatrix{Int}, unit_cell::ExpansionUnitCell)
        coord_index = generate_coord_index(coordinates)
        adj_matrix = zeros(Int, size(coordinates, 2), size(coordinates, 2))
        for (index, col) in enumerate(eachcol(coordinates))
                for bond in unit_cell.bonds
                        if bond.site1 == col[end-1:end]
                                neighbor_coord = neighbor_site(bond, col)
                                ni = get(coord_index, neighbor_coord, nothing)
                                if ni !== nothing
                                        adj_matrix[index, ni] = bond.bond_type
                                        adj_matrix[ni, index] = bond.bond_type
                                end
                        end
                end
        end

        adj_matrix
end

function generate_neighbor_list(coordinates::AbstractMatrix{Int}, unit_cell::UnitCell)
        coord_index = generate_coord_index(coordinates)
        neighbor_list = fill(LatticeVertices{Int}(), size(coordinates, 2))
        for (index, col) in enumerate(eachcol(coordinates))
                for bond in unit_cell.bonds
                        if bond.site1 == col[end]
                                neighbor_index = get(coord_index, neighbor_site(bond, col), nothing)
                                if neighbor_index !== nothing
                                        neighbor_list[index] = union(neighbor_list[index], LatticeVertices(neighbor_index))
                                        neighbor_list[neighbor_index] = union(neighbor_list[neighbor_index], LatticeVertices(index))
                                end
                        end
                end
        end

        neighbor_list
end

function generate_neighbor_list(coordinates::AbstractMatrix{Int}, unit_cell::ExpansionUnitCell)
        coord_index = generate_coord_index(coordinates)
        neighbor_list = fill(ExpansionVertices{Int}(), size(coordinates, 2))
        for (index, col) in enumerate(eachcol(coordinates))
                for bond in unit_cell.expansion_bonds
                        if bond.site1 == col[end]
                                neighbor_index = get(coord_index, neighbor_site(bond, col), nothing)
                                if neighbor_index !== nothing
                                        neighbor_list[index] = union(neighbor_list[index], ExpansionVertices(neighbor_index))
                                        neighbor_list[neighbor_index] = union(neighbor_list[neighbor_index], ExpansionVertices(index))
                                end
                        end
                end
        end

        neighbor_list
end

find_centers(coordinates::AbstractMatrix{Int}) = findall(col -> all(==(0), col[1:end-1]), eachcol(coordinates))

function generate_strong_connections(expansion_coordinates::AbstractMatrix{Int}, lattice_coordinates::AbstractMatrix{Int})

        connections_vec = fill(LatticeVertices{Int}(), size(expansion_coordinates, 2))
        for (i, coord) in enumerate(eachcol(expansion_coordinates))
                connection = findall(==(coord[1:end-1]), eachcol(lattice_coordinates[1:end-2, :]))
                connections_vec[i] = LatticeVertices(collect(connection))
        end
        connections_vec
end








