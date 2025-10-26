struct ClusterExpansionLattice <: AbstractClusterExpansionLattice
    max_order::UInt8
end

centers(lattice::AbstractClusterExpansionLattice) = _NI("centers") # Note that centers should be a subtype of AbstractClusterSet
max_order(lattice::AbstractClusterExpansionLattice) = _NI("max_order")
neighbors(lattice::AbstractClusterExpansionLattice, vs::ExpansionVertices) = _NI("neighbors")
