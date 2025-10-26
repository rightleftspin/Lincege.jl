module GraphHashes

import LINCEGE:
    _NI,
    Vertices.AbstractVertices,
    Lattices.AbstractLattice

abstract type AbstractGraphHash{H<:Unsigned} end

Base.hash(ghash::AbstractGraphHash, h::UInt) = hash(ghash.hash, h)
Base.:(==)(ghash1::AbstractGraphHash, ghash2::AbstractGraphHash) = hash(ghash1) == hash(ghash2)

include("VertexHashes.jl")
include("IsomorphicHashes.jl")
include("TranslationHashes.jl")
include("Permutations.jl")

export AbstractGraphHash,
    VertexHash,
    IsomorphicHash,
    TranslationHash,
    AbstractPermutation,
    IsomorphicPermutation,
    EmptyPermutation,
    ghash,
    mapping
end
