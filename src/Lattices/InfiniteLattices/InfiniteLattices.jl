abstract type AbstractInfiniteLattice <: AbstractLattice end
abstract type AbstractClusterExpansionLattice <: AbstractInfiniteLattice end

n_unique_sites(lattice::AbstractInfiniteLattice) = _NI("n_unique_sites")

n_labels(lattice::AbstractClusterExpansionLattice) = _NI("n_labels")
connections(lattice::AbstractClusterExpansionLattice) = _NI("n_labels")

include("util.jl")
include("SiteExpansionLattices.jl")
include("StrongClusterExpansionLattices.jl")
