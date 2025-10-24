module LCEs

import LINCEGE:
    Vertices.AbstractVertices,
    Lattices.AbstractLattice,
    Clusters.AbstractCluster,
    ClusterSets.AbstractClusterSet,
    _NI

abstract type AbstractLCE{C<:AbstractCluster} end
abstract type AbstractSiteLCE{C<:AbstractCluster} <: AbstractLCE{C} end
abstract type AbstractClusterLCE{C<:AbstractCluster} <: AbstractLCE{C} end

Base.getindex(lce::AbstractLCE, cluster::AbstractCluster, order::Int) = _NI("getindex")
Base.setindex!(lce::AbstractLCE, value, cluster::AbstractCluster, order::Int) = _NI("setindex!")
min_order(lce::AbstractLCE) = _NI("min_order")
add_subgraphs!(lce::AbstractLCE, cluster::AbstractCluster, lattice::AbstractLattice) = _NI("add_subgraphs!")
subgraphs(lce::AbstractLCE, cluster::AbstractCluster) = _NI("subgraphs")

function add_lattice_constant!(lce::AbstractLCE, cluster::AbstractCluster, lattice::AbstractSiteExpansionLattice)
    lce[cluster, length(cluster)] = lattice_constant(cluster)
    add_subgraphs!(lce, cluster, lattice)
end

function add_lattice_constant!(lce::AbstractLCE, cluster::AbstractCluster, lattice::AbstractClusterExpansionLattice)
    lce[cluster, length(cluster) + 1] = lattice_constant(cluster)
    add_subgraphs!(lce, cluster, lattice)
end

function subtract_subgraph_contribution!(lce::AbstractSiteLCE, lattice_constant::Integer, cluster::AbstractCluster, subgraph::AbstractCluster, depth::Integer)
    order = length(cluster)
    lce[subgraph, order] += ((-1) ^ depth) * lattice_constant
    for lower_subgraph in subgraphs(lce, subgraph)
        subtract_subgraph_contribution!(lce, lattice_constant, subgraph, lower_subgraph, depth + 1)
    end
end

function subtract_subgraph_contribution!(lce::AbstractClusterLCE, lattice_constant::Integer, cluster::AbstractCluster, subgraph::AbstractCluster, depth::Integer)
    order = length(cluster) + 1
    lce[subgraph, order] += ((-1) ^ depth) * lattice_constant
    for lower_subgraph in subgraphs(lce, subgraph)
        subtract_subgraph_contribution!(lce, lattice_constant, subgraph, lower_subgraph, depth + 1)
    end
end

function nlce_summation(cluster_set::AbstractClusterSet, lattice::AbstractLattice, lce_type::Type{<:AbstractLCE})
    lce = lce_type(cluster_set, lattice)

    @info "Starting LCE summation for orders $(min_order(lattice)) through $(max_order(lattice))."

    for order in min_order(lattice):max_order(lattice)
        current_order_clusters = get_order(cluster_set, order)
        for cluster in current_order_clusters
            add_lattice_constant!(lce, cluster, lattice)
            for subgraph in subgraphs(lce, cluster)
                subtract_subgraph_contribution!(lce, lattice_constant(cluster), cluster, subgraph, 1)
            end
        end
    end

    @info "LCE summation complete."
    lce
end

end
