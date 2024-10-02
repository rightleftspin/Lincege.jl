# Need to do
#       grow:
#               - add type specifications for the inputs
#               - document the function itself
#               - convert input type for starting vertices to coordinates instead of Int64
#
#       _grow_from_site:
#               - add type specifications for the inputs
#               - document the function itself
#               - optimize the function, it really needs it
#
#       grow_par:
#               - add type specifications for the inputs
#               - document the function itself
#               - convert input type for starting vertices to coordinates instead of Int64
#               - write this function

"""
This is step one of the pipeline. In this step, the algorithm takes in
an underlying cluster, the starting vertices of that cluster and the 
order that clusters should be generated till. Using this information,
the algorithm recursively generates an array of clusters that are
all subclusters of specified order of the underlying cluster.
"""

"""
Main function in step one of the pipeline. Grows clusters from the given
vertices of the underlying_cluster.

Inputs: 
      underlying_cluster: Hashmap connecting vertices to an array of 
      other vertices they are connected to
      
      max_order: Integer that details the order that
      the subclusters will be generated at

      starting_vertices: array of vertex labels (coordinates)
      that are the coordinates that the clusters will be grown from

Output:
      Array of subclusters of the input underlying cluster
"""
function grow(
    underlying_cluster::Union{AbstractNLCECluster,AbstractNLCELattice},
    max_order::Int,
    starting_vertices::Vector{V},
) where {V<:Integer}

    out_array::Vector{AbstractNLCECluster} = Vector()
    guarding_set::Set{V} = Set([])

    for vertex in starting_vertices
        init_neighbors::Set{V} = Set(
            collect(
                filter(
                    neighbor -> !(neighbor in guarding_set),
                    neighbors(underlying_cluster, vertex)[1],
                ),
            ),
        )
        vertices = [vertex]
        _grow_from_site(
            underlying_cluster,
            max_order,
            vertices,
            init_neighbors,
            guarding_set,
            out_array
        )
        push!(guarding_set, vertex)
    end

    out_array
end

"""
Grows the subclusters from a specific site, up till specific order
and outputs them into the out_array.

Inputs: 
      underlying_cluster: Graph with coordinates as vertex labels,
      vertex colors, and edge weights

      max_order: Integer that details the order that
      the subclusters will be generated at

      subclusters_vertices: array of vertices that are the current
      subcluster

      neighbors: set of vertices that are neighbors to the current
      subcluster

      guarding_set: set of vertices not to visit

      out_array: array of subclusters of the underlying_cluster

Output:
      Technically, the output is has_int_leaf, but in practice, the 
      output is the out_array that gets added to.
"""
function _grow_from_site(
    underlying_cluster::Union{AbstractNLCECluster,AbstractNLCELattice},
    max_order::Int,
    subcluster_vertices::Vector{V},
    current_neighbors::Set{V},
    guarding_set::Set{V},
    out_array::Vector{<:AbstractNLCECluster},
) where {V<:Integer}

    if length(subcluster_vertices) == max_order
        push!(out_array, cluster(underlying_cluster, subcluster_vertices))
        return true
    end

    has_int_leaf = false
    new_guarding_set = copy(guarding_set)

    while !isempty(current_neighbors)
        neighbor = pop!(current_neighbors)
        append!(subcluster_vertices, neighbor)

        new_neighbors = copy(current_neighbors)

        for vertex in neighbors(underlying_cluster, neighbor)[1]
            if (
                !(vertex in subcluster_vertices) &
                !(vertex in new_guarding_set) &
                !(vertex in new_neighbors)
            )

                push!(new_neighbors, vertex)
            end
        end

        if _grow_from_site(
            underlying_cluster,
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
        if (nv(underlying_cluster) - length(new_guarding_set)) < max_order
            return (has_int_leaf)
        end
    end
    return (has_int_leaf)
end

"""
Parallel version of the main function in step one of the pipeline.
Grows clusters from the given vertices of the underlying_cluster.

Inputs: 
      underlying_cluster: Graph with coordinates as vertex labels,
      vertex colors, and edge weights

      max_order: Integer that details the order that
      the subclusters will be generated at

      starting_vertices: array of vertex labels (coordinates)
      that are the coordinates that the clusters will be grown from

Output:
      Array of subclusters of the input underlying cluster
"""
function grow_par(underlying_cluster, max_order::Int64, starting_vertices::Vector{Int64}) end
