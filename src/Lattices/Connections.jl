abstract type AbstractConnections end
Base.getindex(c::AbstractConnections, evs::ExpansionVertices{Int}) = _NI("getindex")

struct StrongClusterConnections <: AbstractConnections
        connections::Vector{LatticeVertices{Int}}
end

Base.getindex(c::StrongClusterConnections, evs::ExpansionVertices{Int}) = union(LatticeVertices{Int}(), c.connections[evs])
