struct DirectionCluster{V<:AbstractSet{<:Unsigned},H<:Unsigned} <: AbstractCluster
        expansion_vertices::V
        neighbors::V
        ghash::H
end

function DirectionCluster(vertex::T, lattice::ExpansionLattice, g::DirectionGraph) where {T<:Unsigned}
        DirectionCluster(BitSet(vertex),
                neighbors(lattice, vertex),
                translational_hash(g, real_space_vertices(lattice, vertex), get_mask(lattice, vertex)))
end

function DirectionCluster(vertices::AbstractSet, neighbors::AbstractSet, lattice::ExpansionLattice, g::DirectionGraph)
        DirectionCluster(vertices,
                neighbors,
                translational_hash(g, real_space_vertices(lattice, vertices), get_mask(lattice, vertices)))
end

function neighbor_clusters(cluster::DirectionCluster, lattice::ExpansionLattice, g::DirectionGraph)
        ns = DirectionCluster[]
        for n in cluster.neighbors
                push!(ns, neighbor_cluster(cluster, n, lattice, g))
        end
        ns
end

function neighbor_cluster(cluster::DirectionCluster, n::Unsigned, lattice::ExpansionLattice, g::DirectionGraph)
        DirectionCluster(union(cluster.expansion_vertices, n), union(cluster.neighbors, neighbors(lattice, n)), lattice, g)
end
