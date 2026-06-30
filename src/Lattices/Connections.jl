abstract type AbstractConnections end
Base.getindex(c::AbstractConnections, expansion_vertices::ExpansionVertices{Int}) = _NI("getindex")

struct StrongClusterConnections <: AbstractConnections
        connections::Vector{LatticeVertices{Int}}
end

Base.getindex(c::StrongClusterConnections, expansion_vertices::ExpansionVertices{Int}) = union(LatticeVertices{Int}(), c.connections[expansion_vertices])

struct WeakClusterConnections <: AbstractConnections
        connections::Vector{LatticeVertices{Int}}
        rev_connections::Vector{ExpansionVertices{Int}}
        connections_matrix::Matrix{Int}
        basis_sizes::Vector{Int}
end

function Base.getindex(c::WeakClusterConnections, expansion_vertices::ExpansionVertices{Int})
        lattice_vertices = union(LatticeVertices{Int}(), c.connections[expansion_vertices])
        masking_matrix = .!in.(c.connections_matrix[lattice_vertices, lattice_vertices], (expansion_vertices,))
        masking_matrix[diagind(masking_matrix)] .= 0
        return lattice_vertices, masking_matrix
end

just_lattice_vertices(c::WeakClusterConnections, expansion_vertices::ExpansionVertices{Int}) = union(LatticeVertices{Int}(), c.connections[expansion_vertices])
