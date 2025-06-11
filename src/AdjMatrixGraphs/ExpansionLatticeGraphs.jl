
struct ExpansionLatticeGraph{M<:AbstractMatrix{<:Integer}} <: AbstractAdjMatrixGraph
        adj_matrix::M
end

function ExpansionLatticeGraph(real_space_lattice::RealSpaceLattice, tiling::Tiling)

        dist_matrix = pairwise_distance(real_space_lattice)
        adj_matrix = zeros(Int, size(dist_matrix))

        for ind in findall(d -> any(isapprox.(d, real_space_neighbors(tiling))), dist_matrix)

                @inbounds adj_matrix[ind] = matching_exp_vertex(real_space_lattice, ind[1], ind[2])
        end

        ExpansionLatticeGraph(adj_matrix)
end

function get_mask(g::ExpansionLatticeGraph, expansion_vertices::AbstractVector, real_space_vertices::AbstractVector)
        mask = any.(!in(expansion_vertices), adj_matrix(g)[real_space_vertices, real_space_vertices])
        # Set diagonal to 0 so that it doesn't cancel out later
        mask[diagind(mask)] .= 0
        mask
end
