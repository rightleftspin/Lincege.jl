"""
A julia package for computing LINked Cluster Expansions on a GEneral geometry (LINCEGE).
"""
module LINCEGE

using Base.Threads
using LinearAlgebra
using Distances
using StaticArrays
using NautyGraphs
using JSON3

# Basic not-implemented functionality
include("NI.jl")

# All the nessecary code for constructing a Linked Cluster Expansion
include("Vertices/Vertices.jl")
include("Lattices/Lattices.jl")
include("Clusters/Clusters.jl")
include("ClusterSets/ClusterSets.jl")
include("LCEs/LCEs.jl")

# Extra Physics Related Code, generally slow and not needed for basic LCE construction
include("Physics/Physics.jl")

export Vertices,
    Lattices,
    Clusters,
    ClusterSets,
    LCEs,
    Physics
end
