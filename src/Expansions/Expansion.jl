"""
    Expansion(clusters, lattice, max_order)

Arbitrary expansion of clusters in the NLCE sense, contains the weights necessary to perform the NLCE summation.
"""
struct Expansion <: AbstractExpansion
        expansion_clusters::Dict{UInt,ExpansionCluster}
        order_ids::Vector{Vector{UInt}}
        order_offset::Int
end

function _build_site_expansion(clusters::AbstractClusterSet, lattice::AbstractLattice)
        expansion_clusters = Dict{UInt,ExpansionCluster}()
        sizehint!(expansion_clusters, length(clusters))
        order_ids = [Vector{UInt}() for _ in 1:max_order(lattice)]

        cluster_vec = collect(clusters)
        results = Vector{ExpansionCluster}(undef, length(cluster_vec))
        @threads for i in eachindex(cluster_vec)
                results[i] = ExpansionCluster(cluster_vec[i], clusters, lattice)
        end
        for (cluster, ec) in zip(cluster_vec, results)
                push!(order_ids[length(cluster)], cluster.ghash)
                expansion_clusters[cluster.ghash] = ec
        end

        Expansion(expansion_clusters, order_ids, 0)
end

function _build_cluster_expansion(clusters::AbstractClusterSet, lattice::AbstractLattice)
        expansion_clusters = Dict{UInt,ExpansionCluster}()
        sizehint!(expansion_clusters, length(clusters) + n_unique_sites(clusters))
        order_ids = [Vector{UInt}() for _ in 1:(max_order(lattice)+1)]

        lv::Int = 1
        n_single_site_clusters = n_unique_sites(clusters)
        n_total_sites = length(get_labels(lattice))
        while length(order_ids[1]) < n_single_site_clusters && lv <= n_total_sites

                lv_hash = ghash(clusters, LatticeVertices(lv))

                if !haskey(expansion_clusters, lv_hash)
                        expansion_clusters[lv_hash] = ExpansionCluster(lv, lv_hash, n_single_site_clusters)
                        push!(order_ids[1], lv_hash)
                end

                lv += 1
        end

        cluster_vec = collect(clusters)
        results = Vector{ExpansionCluster}(undef, length(cluster_vec))
        @threads for i in eachindex(cluster_vec)
                results[i] = ExpansionCluster(cluster_vec[i], clusters, lattice)
        end
        for (cluster, ec) in zip(cluster_vec, results)
                push!(order_ids[length(cluster)+1], cluster.ghash)
                expansion_clusters[cluster.ghash] = ec
        end

        Expansion(expansion_clusters, order_ids, 1)
end

function Expansion(clusters::AbstractClusterSet, lattice::SiteExpansionLattice)
        _build_site_expansion(clusters, lattice)
end

function Expansion(clusters::AbstractClusterSet, lattice::FiniteLattice)
        _build_site_expansion(clusters, lattice)
end

function Expansion(clusters::AbstractClusterSet, lattice::AbstractClusterExpansionLattice)
        _build_cluster_expansion(clusters, lattice)
end

function Expansion(clusters::AbstractClusterSet, lattice::AbstractFiniteClusterExpansionLattice)
        _build_cluster_expansion(clusters, lattice)
end

Base.getindex(e::Expansion, cluster_hash::UInt) = e.expansion_clusters[cluster_hash]
each_order(e::Expansion, max_order::Int) = @view e.order_ids[1:max_order]
order_offset(e::Expansion) = e.order_offset

"""
    weights(e::Expansion, order::Int) -> Dict{UInt, Float64}

Return the NLCE weights at the given `order` for all clusters.
Call `summation!` first to populate the underlying weights.
"""
function weights(e::Expansion, order::Int)
        result = Dict{UInt,Float64}()
        for ch in e.order_ids[order]
                cluster = e.expansion_clusters[ch]
                cluster_lattice_constant = cluster.lattice_constant
                for (k, v) in cluster.weights
                        result[k] = get(result, k, 0.0) + cluster_lattice_constant * v
                end
        end
        result
end
