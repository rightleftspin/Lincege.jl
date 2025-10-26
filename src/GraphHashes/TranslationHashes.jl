struct TranslationHash{H<:Unsigned} <: AbstractGraphHash{H}
    hash::H
end

function TranslationHash(lattice::AbstractLattice, vs::AbstractVertices)

end
