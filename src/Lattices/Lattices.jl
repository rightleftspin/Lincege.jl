module Lattices

import LINCEGE:
    _NI,
    Vertices.AbstractVertices,
    Vertices.ExpansionVertices,
    Vertices.LatticeVertices

abstract type AbstractLattice end
abstract type AbstractSiteExpansionLattice <: AbstractLattice end
abstract type AbstractClusterExpansionLattice <: AbstractLattice end
abstract type AbstractRandomLattice <: AbstractLattice end

centers(lattice::AbstractLattice) = _NI("centers")
max_order(lattice::AbstractLattice) = _NI("max_order")
neighbors(lattice::AbstractLattice, vs::AbstractVertices) = _NI("neighbors")

include("SiteExpansionLattices.jl")
include("ClusterExpansionLattices.jl")

export AbstractLattice,
    SiteExpansionLattice,
    ClusterExpansionLattice,
    max_order,
    parallel_dfs
end
