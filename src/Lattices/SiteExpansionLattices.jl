struct SiteExpansionLattice <: AbstractSiteExpansionLattice
    max_order::UInt8
end

centers(lattice::AbstractSiteExpansionLattice) = _NI("centers") # Note that centers should be a subtype of AbstractSiteSet
max_order(lattice::AbstractSiteExpansionLattice) = _NI("max_order")
neighbors(lattice::AbstractSiteExpansionLattice, vs::LatticeVertices) = _NI("neighbors")
