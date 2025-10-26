struct IsomorphicHash{H<:Unsigned} <: AbstractGraphHash{H}
    hash::H
end

function IsomorphicHash(lattice::AbstractLattice, vs::AbstractVertices; perm=false)


end
