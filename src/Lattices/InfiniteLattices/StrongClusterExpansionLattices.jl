struct StrongClusterExpansionLattice <: AbstractClusterExpansionLattice
        max_order::UInt8
        expansion_unit_cell::UnitCell
        expansion_coordinates::Matrix{Int}
        lattice_unit_cells::Vector{UnitCell}
        lattice_coordinates::Vector{Matrix{Int}}
        adj_matrix::Matrix{Int}
        neighbor_list::Vector{ExpansionVertices{Int}}
        connections::StrongClusterConnections
        n_unique_sites::Int
        n_labels::Int
end

function StrongClusterExpansionLattice(max_order::Int, expansion_unit_cell::UnitCell, lattice_unit_cells::Vector{UnitCell})
        @assert max_order > 0 "max_order must be a positive integer"
        @assert length(lattice_unit_cells) == basis_size(expansion_unit_cell) "Need one lattice_unit_cell per expansion site type"

        return StrongClusterExpansionLattice(
                UInt8(max_order),
                expansion_unit_cell,
                lattice_unit_cells,
                phys_coords,
                phys_adj_matrix,
                expansion_neighbor_list,
                StrongClusterConnections(connection_list)
        )
end

centers(lattice::StrongClusterExpansionLattice) = find_centers(lattice.expansion_coordinates)
max_order(lattice::StrongClusterExpansionLattice) = lattice.max_order
n_unique_sites(lattice::StrongClusterExpansionLattice) = lattice.n_unique_sites
n_labels(lattice::StrongClusterExpansionLattice) = lattice.n_labels
neighbors(lattice::StrongClusterExpansionLattice, vs::ExpansionVertices) =
        union(ExpansionVertices(), lattice.neighbor_list[vs])

get_coordinates(lattice::StrongClusterExpansionLattice) = hcat([shift_unit_cell(uc, cs) for (uc, cs) in zip(lattice.lattice_unit_cells, lattice.lattice_coordinates)]...)
get_labels(lattice::StrongClusterExpansionLattice) = vcat([cs[end, :] for cs in lattice.lattice_coordinates]...)
get_site_colors(lattice::StrongClusterExpansionLattice) = vcat([uc.site_colors[cs[end, :]] for (uc, cs) in zip(lattice.lattice_unit_cells, lattice.lattice_coordinates)]...)
bond_matrix(lattice::StrongClusterExpansionLattice) = lattice.adj_matrix
connections(lattice::StrongClusterExpansionLattice) = lattice.connections
