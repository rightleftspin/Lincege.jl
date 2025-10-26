module ClusterExpansions

import LINCEGE:
    Vertices.AbstractVertices,
    Lattices.AbstractLattice,
    Lattices.AbstractSiteExpansionLattice,
    Lattices.AbstractClusterExpansionLattice,
    Clusters.AbstractCluster,
    ClusterSets.AbstractClusters,
    _NI

abstract type AbstractClusterExpansion{C<:AbstractCluster} end
abstract type AbstractSiteClusterExpansion{C<:AbstractCluster} <: AbstractClusterExpansion{C} end
abstract type AbstractClusterClusterExpansion{C<:AbstractCluster} <: AbstractClusterExpansion{C} end

Base.getindex(ce::AbstractClusterExpansion, cluster::AbstractCluster, order::Int) = _NI("getindex")
Base.setindex!(ce::AbstractClusterExpansion, value, cluster::AbstractCluster, order::Int) = _NI("setindex!")
add_subgraphs!(ce::AbstractClusterExpansion, cluster::AbstractCluster, lattice::AbstractLattice) = _NI("add_subgraphs!")
subgraphs(ce::AbstractClusterExpansion, cluster::AbstractCluster) = _NI("subgraphs")

#function add_lattice_constant!(ce::AbstractClusterExpansion, cluster::AbstractCluster, lattice::AbstractSiteExpansionLattice)
#    ce[cluster, length(cluster)] = lattice_constant(cluster)
#    add_subgraphs!(ce, cluster, lattice)
#end
#
#function add_lattice_constant!(LCE::AbstractLCE, cluster::AbstractCluster, lattice::AbstractClusterExpansionLattice)
#    LCE[cluster, length(cluster)+1] = lattice_constant(cluster)
#    add_subgraphs!(LCE, cluster, lattice)
#end
#
#function subtract_subgraph_contribution!(LCE::AbstractSiteLCE, lattice_constant::Integer, cluster::AbstractCluster, subgraph::AbstractCluster, depth::Integer)
#    order = length(cluster)
#    LCE[subgraph, order] += ((-1)^depth) * lattice_constant
#    for lower_subgraph in subgraphs(LCE, subgraph)
#        subtract_subgraph_contribution!(LCE, lattice_constant, subgraph, lower_subgraph, depth + 1)
#    end
#end
#
#function subtract_subgraph_contribution!(LCE::AbstractClusterLCE, lattice_constant::Integer, cluster::AbstractCluster, subgraph::AbstractCluster, depth::Integer)
#    order = length(cluster) + 1
#    LCE[subgraph, order] += ((-1)^depth) * lattice_constant
#    for lower_subgraph in subgraphs(LCE, subgraph)
#        subtract_subgraph_contribution!(LCE, lattice_constant, subgraph, lower_subgraph, depth + 1)
#    end
#end
#
#function nLCE_summation(cluster_set::AbstractClusterSet, lattice::AbstractLattice, LCE_type::Type{<:AbstractLCE})
#    LCE = LCE_type(cluster_set, lattice)
#
#    @info "Starting LCE summation for orders $(min_order(lattice)) through $(max_order(lattice))."
#
#    for order in min_order(lattice):max_order(lattice)
#        current_order_clusters = get_order(cluster_set, order)
#        for cluster in current_order_clusters
#            add_lattice_constant!(LCE, cluster, lattice)
#            for subgraph in subgraphs(LCE, cluster)
#                subtract_subgraph_contribution!(LCE, lattice_constant(cluster), cluster, subgraph, 1)
#            end
#        end
#    end
#
#    @info "LCE summation complete."
#    LCE
#end

end
