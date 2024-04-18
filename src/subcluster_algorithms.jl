"""
    
    subcluster_alg!(cluster::Cluster, max_order::Int, subcluster_vertices::Vector{Int}, neighbors::Set{Int}, guarding_set::Set{Int})

This is a subroutine of the enumerate_subclusters method, it recursively
breaks a cluster into subclusters up till the order specified, sending back everything
by appending to the cluster structure provided.

** Not fully optimized **
"""
function subcluster_alg!(
        cluster::Cluster, 
        max_order::Int, 
        subcluster_vertices::Vector{Int}, 
        neighbors::Set{Int}, 
        guarding_set::Set{Int}
    )
    
    if size(subcluster_vertices)[1] == max_order
        add_subcluster!(cluster, subcluster_vertices)
        return true
    end

    has_int_leaf = false
    new_guarding_set = copy(guarding_set)

    while !isempty(neighbors)
        neighbor = pop!(neighbors)
        append!(subcluster_vertices, neighbor) 

        new_neighbors = copy(neighbors)

        for vertex in outneighbors(cluster, neighbor)
            if (!(vertex in subcluster_vertices) & 
                !(vertex in new_guarding_set) & 
                !(vertex in new_neighbors))

                push!(new_neighbors, vertex)
            end
        end

        if subcluster_alg!(cluster, 
                          max_order, 
                          subcluster_vertices, 
                          new_neighbors, 
                          new_guarding_set)
            pop!(subcluster_vertices)
            has_int_leaf = true
        else
            pop!(subcluster_vertices)
            return(has_int_leaf)
        end
        push!(new_guarding_set, neighbor)
        if (nv(graph(cluster)) - length(new_guarding_set)[1]) < max_order
            return(has_int_leaf)
        end
    end
    return(has_int_leaf)
end

"""

    enumerate_subclusters!(cluster::Cluster, max_order::Int, starting_vertices::Vector{Int})

Takes in a specific cluster and adds to the cluster, every subcluster of it
up until the provided order and starting from the starting vertices provided
"""
function enumerate_subclusters!(cluster::Cluster, max_order::Int, starting_vertices::Vector{Int})
    
    guarding_set::Set{Int} = Set([])

    for vertex in starting_vertices
        neighbors::Set{Int} = Set(collect(filter(neighbor -> !(neighbor in guarding_set), collect(outneighbors(cluster, vertex)))))
        vertices = [vertex]
        subcluster_alg!(
            cluster, 
            max_order, 
            vertices, 
            neighbors, 
            guarding_set
        )
        push!(guarding_set, vertex)
    end
end

"""
    
    subcluster_alg!(cluster::Cluster, max_order::Int, subcluster_vertices::Vector{Int}, neighbors::Set{Int}, guarding_set::Set{Int})

This is a subroutine of the enumerate_subclusters method, it recursively
breaks a cluster into subclusters up till the order specified, sending back everything
by appending to the cluster structure provided.

** Not fully optimized **
"""
function subcluster_alg_full!(
        cluster::Cluster, 
        max_order::Int, 
        subcluster_vertices::Vector{Int}, 
        neighbors::Set{Int}, 
        guarding_set::Set{Int}
    )

    add_subcluster!(cluster, subcluster_vertices)

    if size(subcluster_vertices)[1] == max_order
        return true
    end

    has_int_leaf = false
    new_guarding_set = copy(guarding_set)

    while !isempty(neighbors)
        neighbor = pop!(neighbors)
        append!(subcluster_vertices, neighbor) 

        new_neighbors = copy(neighbors)

        for vertex in outneighbors(cluster, neighbor)
            if (!(vertex in subcluster_vertices) & 
                !(vertex in new_guarding_set) & 
                !(vertex in new_neighbors))

                push!(new_neighbors, vertex)
            end
        end

        if subcluster_alg_full!(cluster, 
                          max_order, 
                          subcluster_vertices, 
                          new_neighbors, 
                          new_guarding_set)
            pop!(subcluster_vertices)
            has_int_leaf = true
        else
            pop!(subcluster_vertices)
            return(has_int_leaf)
        end
        push!(new_guarding_set, neighbor)
        if (nv(graph(cluster)) - length(new_guarding_set)[1]) < max_order
            return(has_int_leaf)
        end
    end
    return(has_int_leaf)
end

"""

    enumerate_subclusters!(cluster::Cluster, max_order::Int, starting_vertices::Vector{Int})

Takes in a specific cluster and adds to the cluster, every subcluster of it
up until the provided order and starting from the starting vertices provided
"""
function enumerate_subclusters_full!(cluster::Cluster, max_order::Int, starting_vertices::Vector{Int})
    
    guarding_set::Set{Int} = Set([])

    for vertex in starting_vertices
        neighbors::Set{Int} = Set(collect(filter(neighbor -> !(neighbor in guarding_set), collect(outneighbors(cluster, vertex)))))
        starting_subcluster_vertices = [vertex]
        subcluster_alg_full!(
            cluster, 
            max_order, 
            starting_subcluster_vertices, 
            neighbors, 
            guarding_set
        )
        push!(guarding_set, vertex)
    end
end
