
abstract type AbstractBundle end


"""
This is the cluster struct that powers the entire NLCE algorithm. It is designed to be very general so
it deals with many use cases.
"""
struct SiteExpansionBundle <: AbstractBundle

    "Underlying lattice for this NLCE expansion"
    lattice::AbstractCluster
    "Set full of all translationally invariant clusters"
    clusters::AbstractSet{<:AbstractCluster}


    function SiteExpansionBundle(
    )

        return new{has_underlying_cluster, cluster_expansion, vertex_labeled, edge_labeled}(
            coordinates,
            start,
            adj_list,
            adj_matrices,
            underlying_cluster,
            vertices,
            coordinate_bundles,
            super_adj_list,
        )
    end
end
