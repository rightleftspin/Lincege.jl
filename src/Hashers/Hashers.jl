module Hashers

using Base.Threads
using LinearAlgebra
using NautyGraphs

import LINCEGE:
        _NI,
        Vertices.ExpansionVertices,
        Vertices.LatticeVertices,
        Lattices.get_coordinates,
        Lattices.get_labels,
        Lattices.get_site_colors,
        Lattices.bond_matrix,
        Lattices.AbstractInfiniteLattice,
        Lattices.AbstractClusterExpansionLattice,
        Lattices.AbstractConnections

abstract type AbstractHasher end

ghash(h::AbstractHasher, evs::ExpansionVertices) = _NI("ghash")
ghash(h::AbstractHasher, lvs::LatticeVertices) = _NI("ghash")

include("util.jl")
include("TranslationHasher.jl")
include("IsomorphicHasher.jl")
include("SymmetricHasher.jl")

end
