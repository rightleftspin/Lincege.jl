"""
    Expansion(clusters, lattice, max_order)

Arbitrary expansion of clusters in the NLCE sense, contains the weights necessary to perform the NLCE summation.
"""
struct Expansion <: AbstractExpansion
        expansion_clusters::Dict{UInt,ExpansionCluster}
        order_ids::Vector{Vector{UInt}}
        order_offset::Int
end

function Expansion(clusters::AbstractClusterSet, lattice::SiteExpansionLattice)
        expansion_clusters = Dict{UInt,ExpansionCluster}()
        sizehint!(expansion_clusters, length(clusters))
        order_ids = [Vector{UInt}() for _ in 1:max_order(lattice)]

        for cluster in clusters
                push!(order_ids[length(cluster)], cluster.ghash)
                expansion_clusters[cluster.ghash] = ExpansionCluster(cluster, clusters, lattice)
        end

        Expansion(expansion_clusters, order_ids, 0)
end

function Expansion(clusters::AbstractClusterSet, lattice::AbstractClusterExpansionLattice)
        expansion_clusters = Dict{UInt,ExpansionCluster}()
        sizehint!(expansion_clusters, length(clusters) + n_unique_sites(clusters))
        order_ids = [Vector{UInt}() for _ in 1:(max_order(lattice)+1)]

        # Adds Single Sites to the cluster expansion
        lv::Int = 1
        n_single_site_clusters = n_unique_sites(clusters)
        while length(order_ids[1]) < n_single_site_clusters

                lv_hash = ghash(clusters, LatticeVertices(lv))

                if !haskey(expansion_clusters, lv_hash)
                        expansion_clusters[lv_hash] = ExpansionCluster(lv, lv_hash, n_single_site_clusters)
                        push!(order_ids[1], lv_hash)
                end

                lv += 1
        end

        for cluster in clusters
                push!(order_ids[length(cluster)+1], cluster.ghash)
                expansion_clusters[cluster.ghash] = ExpansionCluster(cluster, clusters, lattice)
        end

        Expansion(expansion_clusters, order_ids, 1)
end

Base.getindex(e::Expansion, cluster_hash::UInt) = e.expansion_clusters[cluster_hash]
each_order(e::Expansion, max_order::Int) = e.order_ids[1:max_order]
order_offset(e::Expansion) = e.order_offset

function weights(e::Expansion, order::Int)
        result = Dict{UInt,Float64}()
        for ch in e.order_ids[order]
                cluster = e.expansion_clusters[ch]
                lc = cluster.lattice_constant
                for (k, v) in cluster.weights
                        result[k] = get(result, k, 0.0) + lc * v
                end
        end
        result
end

"""
    write_to_json(expansion, lattice, filepath)

Serialize an `Expansion` and its associated `lattice` geometry to a JSON file at
`filepath`.  The file contains a JSON array — one object per cluster — with the
following fields:

- `cluster_hash`   — unique cluster identifier (string representation of UInt64)
- `order`          — 1-based order index into the expansion's `order_ids`
- `n_sites`        — number of sites in the cluster
- `coordinates`    — list of Cartesian coordinate vectors (one per site)
- `site_colors`    — list of integer site-color labels (one per site)
- `bonds`          — edge list `[i, j, weight]` using local 1-based indices
- `weights`        — vector of weights per order for the cluster
"""
function write_to_json(e::Expansion, lattice::AbstractLattice, filepath::String)
        all_coords = get_coordinates(lattice)
        all_colors = get_site_colors(lattice)
        adj = bond_matrix(lattice)
        all_weights = [weights(e, i) for i in 1:length(e.order_ids)]

        clusters_data = Vector{Dict{String,Any}}()

        for (order_idx, cluster_hashes) in enumerate(e.order_ids)
                for cluster_hash in cluster_hashes
                        cluster = e.expansion_clusters[cluster_hash]
                        vs = cluster.vertices
                        n = length(vs)

                        coords = [collect(col) for col in eachcol(all_coords[:, vs])]
                        colors = collect(all_colors[vs])
                        bonds = [[b[1], b[2], b[3]] for b in adj_mat_to_edge_list(adj[vs, vs])]
                        wts = [get(d, cluster_hash, 0.0) for d in all_weights]

                        push!(clusters_data, Dict(
                                "cluster_hash" => string(cluster_hash),
                                "order" => order_idx,
                                "n_sites" => n,
                                "coordinates" => coords,
                                "site_colors" => colors,
                                "bonds" => bonds,
                                "weights" => wts
                        ))
                end
        end

        open(filepath, "w") do io
                JSON.print(io, clusters_data, 2)
        end
end

function _expansion_table_data(e::Expansion, cs::Vector{<:AbstractClusterSet}, max_order::Int)
        off = order_offset(e)
        table_rows = Vector{Vector{Any}}()
        all_weights = [weights(e, i) for i in 1:length(e.order_ids)]

        for order in (1-off):max_order
                idx = order + off
                row = Vector{Any}()
                push!(row, order)
                if order != 0
                        for c in cs
                                push!(row, count(x -> length(x) == order, c))
                        end
                else
                        for _ in cs
                                push!(row, 1)
                        end

                end

                lc_sum = 0
                sg_sum = 0
                for id in e.order_ids[idx]
                        lc_sum += e.expansion_clusters[id].lattice_constant
                        sg_sum += length(e.expansion_clusters[id].subgraphs)
                end
                push!(row, lc_sum)
                push!(row, sg_sum)
                push!(row, sum(p -> p.second, all_weights[idx]))

                push!(table_rows, row)
        end
        data = permutedims(reduce(hcat, table_rows))
        cs_labels = [_cs_column_label(c) for c in cs]
        return data, cs_labels
end

