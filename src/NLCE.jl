module NLCE

using NautyGraphs
using StaticArrays

# Add the relevant structs
include("Clusters.jl")

# Add the relevant helper functions
include("helpers/util.jl")
include("helpers/wrappers.jl")
include("helpers/writers.jl")
include("helpers/symmetries.jl")

# Add the basic pipeline
include("pipeline/prune.jl")
include("pipeline/propagate.jl")
include("pipeline/combine.jl")

export simple_NLCE, coord_NLCE, write_to_file, write_to_file_fortran

end
