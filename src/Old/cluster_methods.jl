function reindex_to_subcluster(
        cluster::Cluster,
        super_vertices::AbstractVector{<:Integer},
)
        # Sort them for translational hashing
        sorted_vertices = sort(unique(vcat(connections(cluster)[super_vertices]...)))

        new_adj_list = reindex_adj_list(adj_list(cluster), super_vertices)
        new_connections =
                reindex_connections(connections(cluster), super_vertices, sorted_vertices)
        new_adjacency_matrices = reindex_adjacency_matrices(
                adjacency_matrices(cluster),
                super_vertices,
                sorted_vertices,
        )

        new_adj_list, new_connections, new_adjacency_matrices
end

