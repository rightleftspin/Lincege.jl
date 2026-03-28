function get_subgraphs(c::AbstractCluster, lattice::AbstractLattice)
        if length(c) == 1
                return Set()
        end

        max_depth = length(c) - 1
        ctrs = c.vs
        roots = [typeof(ctrs)(center) for center in ctrs]
        visited = Set()

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
                        if v in c.vs
                                dfs(union(cluster, typeof(ctrs)(v)))
                        end
                end
        end

        for c in roots
                dfs(c)
        end

        visited
end

function adj_mat_to_edge_list(adj_matrix::AbstractMatrix{<:Real})
        edge_list = []

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
