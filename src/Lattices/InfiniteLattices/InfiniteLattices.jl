abstract type AbstractInfiniteLattice <: AbstractLattice end

num_unique_sites(lattice::AbstractInfiniteLattice) = _NI("n_unique_sites")

include("util.jl")
include("SiteExpansionLattices.jl")
