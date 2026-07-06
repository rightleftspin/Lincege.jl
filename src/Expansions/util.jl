function get_subgraphs(c::AbstractCluster, lattice::AbstractLattice)
        cluster_vertices = c.vertices
        if length(c) == 1
                return Set{typeof(cluster_vertices)}()
        end

        max_depth = length(c) - 1
        roots = [typeof(cluster_vertices)(center) for center in cluster_vertices]
        visited = Set{typeof(cluster_vertices)}()

        function try_mark(cluster)
                already = cluster in visited
                if !already
                        push!(visited, cluster)
                end

                if length(cluster) == max_depth
                        return false
                end
                !already
        end

        function dfs(cluster)
                if !try_mark(cluster)
                        return
                end
                for v in neighbors(lattice, cluster)
                        if v in c.vertices
                                dfs(union(cluster, typeof(cluster_vertices)(v)))
                        end
                end
        end

        for root in roots
                dfs(root)
        end

        visited
end

function adj_mat_to_edge_list(adj_matrix::AbstractMatrix{<:Real})
        edge_list = Tuple{Int,Int,Int}[]

        n = size(adj_matrix, 1)
        for i in 1:n
                for j in i+1:n
                        weight = adj_matrix[i, j]
                        if weight != 0
                                push!(edge_list, (i, j, weight))
                        end
                end
        end

        edge_list
end

