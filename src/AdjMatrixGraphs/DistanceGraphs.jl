struct DistanceGraph{W<:AbstractVector,M<:AbstractMatrix{<:Integer}} <: AbstractAdjMatrixGraph
        weight_info::W
        adj_matrix::M
end

function DistanceGraph(real_space_lattice::RealSpaceLattice)

        dist_matrix = pairwise_distance(real_space_lattice)
        adj_matrix = zeros(Int, size(dist_matrix))
        neighbors = unique(dist_matrix)

        for (i, neighbor) in enumerate(neighbors)
                @inbounds adj_matrix[findall(isapprox(neighbor), dist_matrix)] .= 2^i
        end

        DistanceGraph(neighbors, adj_matrix)
end

distance_hash(g::DistanceGraph, real_space_vertices::AbstractVector) = hash(sum(g.adj_matrix[real_space_vertices, real_space_vertices], dims=2))
