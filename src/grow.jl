using Graphs
using DataStructures
using StaticArrays

begin # Growing functions to get subclusters from a cluster

        struct Subgraph{T<:Integer}
                vertices::AbstractSet{T} # needs to be a sorted container since it is in order
                neighbors::AbstractSet{T}
        end

        function Subgraph(vertex::T, cluster::Cluster) where {T<:Integer}
                Subgraph{T}(BitSet([vertex]), BitSet(neighbors(cluster, vertex)))
        end

        function neighbor_subgraphs(subgraph::Subgraph, cluster::Cluster)
                neighbors = Vector{Subgraph}(undef, length(subgraph.neighbors))
                neighbor_subgraphs!(subgraph, cluster, neighbors)
        end

        function neighbor_subgraphs!(subgraph, cluster, ns)
                for (i, neighbor) in enumerate(subgraph.neighbors)
                        new_vertices = union(subgraph.vertices, neighbor)
                        new_neighbors = setdiff(union(subgraph.neighbors, neighbors(cluster, neighbor)), subgraph.vertices)
                        ns[i] = Subgraph(new_vertices, new_neighbors)
                end
                ns
        end

        to_vector(subgraph::Subgraph) = collect(subgraph.vertices)
        nodes(subgraph::Subgraph) = subgraph.vertices

        Base.length(subgraph::Subgraph) = length(subgraph.vertices)

        Base.hash(subgraph::Subgraph, h::UInt) = hash(subgraph.vertices, h)
        Base.isequal(g1::Subgraph, g2::Subgraph) = (g1.vertices == g2.vertices)


        struct Consensus
                is_done::AbstractVector
                thread_id::Int
        end

        all_threads_done(c::Consensus) = (sum(c.is_done) == 0)
        thread_empty!(c::Consensus) = (c.is_done[c.thread_id] = 0)
        thread_working!(c::Consensus) = (c.is_done[c.thread_id] = 1)


        function dfs_kernel!(dfs_lock::ReentrantLock, next::Stack{Subgraph}, cluster::Cluster, parents::AbstractSet{Subgraph}, max_order::Integer, c::Consensus)

                while true

                        lock(dfs_lock)
                        if all_threads_done(c)
                                unlock(dfs_lock)
                                break
                        end

                        if isempty(next)
                                println("thread $(c.thread_id) empty")
                                thread_empty!(c)
                                println("$(c.is_done) state")
                                unlock(dfs_lock)
                                sleep(0.1)
                                continue
                        end

                        thread_working!(c)
                        subgraph = pop!(next)
                        println("thread $(c.thread_id) working")

                        unlock(dfs_lock)
                        if length(subgraph) < max_order
                                filtered_subgraphs = filter(!in(parents), neighbor_subgraphs(subgraph, cluster))
                                lock(dfs_lock)
                                for fs in filtered_subgraphs
                                        push!(parents, fs)
                                        push!(next, fs)
                                end
                                unlock(dfs_lock)

                        end
                end
                println("Finished")
        end

        function dfs!(next::Stack{Subgraph}, cluster::Cluster, source::Subgraph, parents::AbstractSet{Subgraph}, max_order::Integer)

                push!(next, source)
                push!(parents, source)

                dfs_lock = ReentrantLock()

                con_array = ones(Int, Threads.nthreads())

                Threads.@threads for i in 1:Threads.nthreads()
                        c = Consensus(view(con_array, :), Threads.threadid())
                        dfs_kernel!(dfs_lock, next, cluster, parents, max_order, c)
                end

                parents
        end

        function grow_lower_site!(cluster::Cluster, parents::AbstractSet{Subgraph}, start::Integer, max_order::Integer)

                root_subgraph = Subgraph(start, cluster)

                next = Stack{Subgraph}()

                dfs!(next, cluster, root_subgraph, parents, max_order)
        end


        function grow_lower(cluster::Cluster, start::Vector{<:Integer}, max_order::Integer)

                parents = Set{Subgraph}()

                for vertex in start
                        grow_lower_site!(cluster, parents, vertex, max_order)
                end

                println("here")
                to_vector.(parents)
        end

        function grow_lower_rec(cluster::Cluster, start::Vector{<:Integer}, max_order::Integer)

                out_array::Vector{Vector{Int}} = Vector()
                guarding_set::Set{Int} = Set([])

                for vertex in start
                        init_neighbors::Set{Int} = Set(
                                collect(
                                        filter(
                                                neighbor -> !(neighbor in guarding_set),
                                                neighbors(cluster, vertex),
                                        ),
                                ),
                        )
                        vertices = [vertex]
                        _grow_from_site_lower(
                                cluster,
                                max_order,
                                vertices,
                                init_neighbors,
                                guarding_set,
                                out_array,
                        )
                        push!(guarding_set, vertex)
                end

                out_array
        end

        """
        This is step one of the NLCE pipeline. The algorithm takes in
        a lattice and the order that clusters should be generated till.
        Using this information, the algorithm recursively generates an array of
        clusters that are all subclusters of specified order or lower of the lattice.
        """
        function grow_exact(cluster::Cluster, start::Vector{<:Integer}, max_order::Integer)

                out_array::Vector{Vector{Int}} = Vector()
                guarding_set::Set{Int} = Set([])

                for vertex in start
                        init_neighbors::Set{Int} = Set(
                                collect(
                                        filter(
                                                neighbor -> !(neighbor in guarding_set),
                                                neighbors(cluster, vertex),
                                        ),
                                ),
                        )
                        vertices = [vertex]
                        _grow_from_site_exact(
                                cluster,
                                max_order,
                                vertices,
                                init_neighbors,
                                guarding_set,
                                out_array,
                        )
                        push!(guarding_set, vertex)
                end

                out_array
        end

        """
        Grows the subclusters from a specific site, up till specific order and outputs them into
        the out_array. This adds all subclusters at or below the given max_order.
        """
        function _grow_from_site_lower(
                cluster::Cluster,
                max_order::Integer,
                subcluster_vertices::AbstractVector{V},
                current_neighbors::Set{V},
                guarding_set::Set{V},
                out_array::AbstractVector{<:AbstractVector{V}},
        ) where {V<:Integer}

                push!(out_array, deepcopy(subcluster_vertices))

                if length(subcluster_vertices) == max_order
                        return true
                end

                has_int_leaf = false
                new_guarding_set = copy(guarding_set)

                while !isempty(current_neighbors)
                        neighbor = pop!(current_neighbors)
                        append!(subcluster_vertices, neighbor)

                        new_neighbors = copy(current_neighbors)

                        for vertex in neighbors(cluster, neighbor)
                                if (
                                        !(vertex in subcluster_vertices) & !(vertex in new_guarding_set) &
                                        !(vertex in new_neighbors)
                                )

                                        push!(new_neighbors, vertex)
                                end
                        end

                        if _grow_from_site_lower(
                                cluster,
                                max_order,
                                subcluster_vertices,
                                new_neighbors,
                                new_guarding_set,
                                out_array,
                        )
                                pop!(subcluster_vertices)
                                has_int_leaf = true
                        else
                                pop!(subcluster_vertices)
                                return (has_int_leaf)
                        end
                        push!(new_guarding_set, neighbor)

                        if (nsv(cluster) - length(new_guarding_set)) < max_order
                                return (has_int_leaf)
                        end
                end
                return (has_int_leaf)
        end

        """
        Grows the subclusters from a specific site, up till specific order and outputs them into
        the out_array. This adds all subclusters at the given max_order.
        """
        function _grow_from_site_exact(
                cluster::Cluster,
                max_order::Integer,
                subcluster_vertices::AbstractVector{V},
                current_neighbors::Set{V},
                guarding_set::Set{V},
                out_array::AbstractVector{<:AbstractVector{V}},
        ) where {V<:Integer}


                if length(subcluster_vertices) == max_order
                        push!(out_array, deepcopy(subcluster_vertices))
                        return true
                end

                has_int_leaf = false
                new_guarding_set = copy(guarding_set)

                while !isempty(current_neighbors)
                        neighbor = pop!(current_neighbors)
                        append!(subcluster_vertices, neighbor)

                        new_neighbors = copy(current_neighbors)

                        for vertex in neighbors(cluster, neighbor)
                                if (
                                        !(vertex in subcluster_vertices) & !(vertex in new_guarding_set) &
                                        !(vertex in new_neighbors)
                                )

                                        push!(new_neighbors, vertex)
                                end
                        end

                        if _grow_from_site_exact(
                                cluster,
                                max_order,
                                subcluster_vertices,
                                new_neighbors,
                                new_guarding_set,
                                out_array,
                        )
                                pop!(subcluster_vertices)
                                has_int_leaf = true
                        else
                                pop!(subcluster_vertices)
                                return (has_int_leaf)
                        end
                        push!(new_guarding_set, neighbor)

                        if (nsv(cluster) - length(new_guarding_set)) < max_order
                                return (has_int_leaf)
                        end
                end
                return (has_int_leaf)
        end
end
