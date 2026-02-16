struct ClusterSet{C,H}
    clusters::Set{C}
    hasher::H
end

function ClusterSet(lattice::AbstractLattice)
    ClusterSet{Cluster,TranslationHasher}(
        Set{Cluster}(),
        TranslationHasher(lattice)
    )
end
