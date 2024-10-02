# Need to do
#       NLCELattice:
#                       - Raise errors if adj_list and adj_matrix don't make sense
"""
"""

abstract type AbstractNLCECluster end

struct NLCECluster <: AbstractNLCECluster
    vertex_labels::SubArray
    adj_list::AbstractVector{<:SubArray}
    adj_matrix::SubArray
    adj_list_weights::AbstractVector{<:AbstractVector{<:SubArray}}
    adj_matrix_weights::AbstractVector{<:SubArray}
    
    # Basic Constructor
    function NLCECluster(
        vertex_labels::SubArray,
        adj_list::AbstractVector{<:SubArray},
        adj_matrix::SubArray,
        adj_list_weights::AbstractVector{<:AbstractVector{<:SubArray}},
        adj_matrix_weights::AbstractVector{<:SubArray}
    )

        return new(
            vertex_labels,
            adj_list,
            adj_matrix,
            adj_list_weights,
            adj_matrix_weights,
        )
    end
end

begin #Required functions for the pipeline
    nv(cluster::NLCECluster) = length(cluster.vertex_labels)
    neighbors(cluster::NLCECluster, vertex::Int64) = @view cluster.adj_list[vertex]
end
