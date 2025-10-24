struct DirectionClusters{C<:DirectionCluster} <: AbstractClusterSet{C}
    clusters::Dict{<:Unsigned, C}
end
