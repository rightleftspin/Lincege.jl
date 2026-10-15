"""
    write_to_json(expansion, lattice, filepath)

Serialize an `Expansion` and its associated `lattice` geometry to a JSON file at
`filepath`. 

- `cluster_hash`   - unique cluster identifier (string representation of UInt64)
- `order`          - 1-based order index into the expansion's `order_ids`
- `n_sites`        - number of sites in the cluster
- `coordinates`    - list of Cartesian coordinate vectors (one per site)
- `site_colors`    - list of integer site-color labels (one per site)
- `bonds`          - edge list `[i, j, weight]` using local 1-based indices
- `weights`        - vector of weights per order for the cluster
"""
function write_to_json(e::Expansion, lattice::AbstractLattice, filepath::String)
        data = _expansion_export_data(e, lattice)
        open(filepath, "w") do io
                JSON.print(io, data, 2)
        end
end

function _expansion_export_data(e::Expansion, lattice::AbstractLattice)
        all_coords = get_coordinates(lattice)
        all_colors = get_site_colors(lattice)
        adj = bond_matrix(lattice)
        all_weights = [weights(e, i) for i in 1:length(e.order_ids)]

        clusters_data = Vector{Dict{String,Any}}()

        for (order_idx, cluster_hashes) in enumerate(e.order_ids)
                for cluster_hash in cluster_hashes
                        cluster = e.expansion_clusters[cluster_hash]
                        vertices = cluster.vertices
                        n = length(vertices)

                        coords = [collect(col) for col in eachcol(all_coords[:, vertices])]
                        colors = collect(all_colors[vertices])
                        bonds = [[b[1], b[2], b[3]] for b in adj_mat_to_edge_list(adj[vertices, vertices])]
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

        return clusters_data
end


"""
    write_to_json(iso_expansion, iso_set, sym_expansion, sym_set,  lattice, filepath)

Serialize an `Expansion` and its associated `lattice` geometry to a JSON file at
`filepath`. The file contains a similar output to the original write_to_json function,
but with the addition of relating multiple symmetrically distinct clusters to their
isomorphic equivalent.
"""
function write_to_json(iso_expansion::Expansion, iso_set::ClusterSet{Ci,IsomorphicHasher},
        sym_expansion::Expansion, sym_set::ClusterSet{Cs,SymmetricHasher},
        lattice::AbstractLattice, filepath::String) where {Ci,Cs}
        data = _expansion_export_data(iso_expansion, iso_set, sym_expansion, sym_set, lattice)
        open(filepath, "w") do io
                JSON.print(io, data, 2)
        end
end

function _expansion_export_data(iso_expansion::Expansion, iso_set::ClusterSet{Ci,IsomorphicHasher},
        sym_expansion::Expansion, sym_set::ClusterSet{Cs,SymmetricHasher},
        lattice::AbstractLattice) where {Ci,Cs}

        iso_hasher = iso_set.hasher
        all_coords = get_coordinates(lattice)
        all_colors = get_site_colors(lattice)
        adj = bond_matrix(lattice)
        iso_weights = [weights(iso_expansion, i) for i in 1:length(iso_expansion.order_ids)]
        sym_weights = [weights(sym_expansion, i) for i in 1:length(sym_expansion.order_ids)]

        entries = Dict{UInt,Dict{String,Any}}()

        for (order_idx, cluster_hashes) in enumerate(iso_expansion.order_ids)
                for cluster_hash in cluster_hashes
                        cluster = iso_expansion.expansion_clusters[cluster_hash]
                        vertices = cluster.vertices
                        n = length(vertices)

                        _, perm = ghash_with_permutation(iso_hasher, vertices)
                        order = sortperm(collect(perm))

                        sub_adj = adj[vertices, vertices]
                        sub_colors = collect(all_colors[vertices])
                        canonical_adj = sub_adj[order, order]
                        canonical_colors = sub_colors[order]

                        canonical_bonds = [[b[1], b[2], b[3]] for b in adj_mat_to_edge_list(canonical_adj)]

                        entries[cluster_hash] = Dict{String,Any}(
                                "cluster_hash" => string(cluster_hash),
                                "order" => order_idx,
                                "n_sites" => n,
                                "canonical_site_colors" => canonical_colors,
                                "canonical_bonds" => canonical_bonds,
                                "weights" => [get(d, cluster_hash, 0.0) for d in iso_weights],
                                "related_symmetric" => Vector{Dict{String,Any}}()
                        )
                end
        end

        for cluster_hashes in sym_expansion.order_ids
                for cluster_hash in cluster_hashes
                        cluster = sym_expansion.expansion_clusters[cluster_hash]
                        vertices = cluster.vertices

                        iso_hash, perm = ghash_with_permutation(iso_hasher, vertices)
                        haskey(entries, iso_hash) || error("Symmetric cluster $(cluster_hash) has no related isomorphic cluster (iso hash $(iso_hash)); both expansions must be built over the same lattice.")

                        order = sortperm(collect(perm))
                        permuted_vertices = collect(vertices)[order]
                        push!(entries[iso_hash]["related_symmetric"], Dict{String,Any}(
                                "cluster_hash" => string(cluster_hash),
                                "permuted_coordinates" => [collect(col) for col in eachcol(all_coords[:, permuted_vertices])],
                                "weights" => [get(d, cluster_hash, 0.0) for d in sym_weights]
                        ))
                end
        end

        return collect(values(entries))
end

_txt_bonds(bonds) = join((join(bond, ' ') for bond in bonds), ' ')
_txt_coordinates(coordinates) = join(("(" * join(coordinate, ',') * ")" for coordinate in coordinates), ' ')
_write_txt_record(io::IO, fields...) = println(io, join(fields, " : "))

"""
    write_to_txt_file(expansion, lattice, filepath)

Write the same cluster data as `write_to_json` in a colon-separated text format,
with one line per cluster:

`n_sites : cluster_hash : bonds : coordinates : site_colors : weights : order`
"""
function write_to_txt_file(e::Expansion, lattice::AbstractLattice, filepath::AbstractString)
        data = _expansion_export_data(e, lattice)
        sort!(data; by=entry -> (entry["order"], parse(UInt, entry["cluster_hash"])))
        open(filepath, "w") do io
                println(io, "# n_sites : cluster_hash : bonds : coordinates : site_colors : weights : order")
                for entry in data
                        _write_txt_record(io,
                                entry["n_sites"], entry["cluster_hash"],
                                _txt_bonds(entry["bonds"]),
                                _txt_coordinates(entry["coordinates"]),
                                join(entry["site_colors"], ' '),
                                join(entry["weights"], ' '), entry["order"],
                        )
                end
        end
end

"""
    write_to_txt_file(iso_expansion, iso_set, sym_expansion, sym_set, lattice, filepath)

Write all fields of the joint `write_to_json` export as colon-separated records:

`isomorphic : n_sites : cluster_hash : canonical_bonds : canonical_site_colors : weights : order`

`symmetric : isomorphic_hash : cluster_hash : permuted_coordinates : weights`
"""
function write_to_txt_file(iso_expansion::Expansion, iso_set::ClusterSet{Ci,IsomorphicHasher},
        sym_expansion::Expansion, sym_set::ClusterSet{Cs,SymmetricHasher},
        lattice::AbstractLattice, filepath::AbstractString) where {Ci,Cs}
        data = _expansion_export_data(iso_expansion, iso_set, sym_expansion, sym_set, lattice)
        sort!(data; by=entry -> (entry["order"], parse(UInt, entry["cluster_hash"])))
        open(filepath, "w") do io
                println(io, "# isomorphic : n_sites : cluster_hash : canonical_bonds : canonical_site_colors : weights : order")
                println(io, "# symmetric : isomorphic_hash : cluster_hash : permuted_coordinates : weights")
                for entry in data
                        _write_txt_record(io,
                                "isomorphic", entry["n_sites"], entry["cluster_hash"],
                                _txt_bonds(entry["canonical_bonds"]),
                                join(entry["canonical_site_colors"], ' '),
                                join(entry["weights"], ' '), entry["order"],
                        )
                        related = sort(entry["related_symmetric"]; by=cluster -> parse(UInt, cluster["cluster_hash"]))
                        for cluster in related
                                _write_txt_record(io,
                                        "symmetric", entry["cluster_hash"], cluster["cluster_hash"],
                                        _txt_coordinates(cluster["permuted_coordinates"]),
                                        join(cluster["weights"], ' '),
                                )
                        end
                end
        end
end

_expansion_table_data(e::AbstractExpansion, cluster_sets::Vector{<:AbstractClusterSet}, max_order::Int) = _NI("_expansion_table_data")

function _expansion_table_data(e::Expansion, cluster_sets::Vector{<:AbstractClusterSet}, max_order::Int)
        off = order_offset(e)
        table_rows = Vector{Vector{Any}}()
        all_weights = [weights(e, i) for i in 1:length(e.order_ids)]

        for order in (1-off):max_order
                idx = order + off
                row = Vector{Any}()
                push!(row, order)
                if order != 0
                        for c in cluster_sets
                                push!(row, count(x -> length(x) == order, c))
                        end
                else
                        for _ in cluster_sets
                                push!(row, 1)
                        end

                end

                lattice_constant_sum = 0
                subgraph_sum = 0
                for id in e.order_ids[idx]
                        lattice_constant_sum += e.expansion_clusters[id].lattice_constant
                        subgraph_sum += length(e.expansion_clusters[id].subgraphs)
                end
                push!(row, lattice_constant_sum)
                push!(row, subgraph_sum)
                push!(row, sum(p -> p.second, all_weights[idx]))

                push!(table_rows, row)
        end
        data = permutedims(reduce(hcat, table_rows))
        cluster_set_labels = [_cluster_set_column_label(c) for c in cluster_sets]
        return data, cluster_set_labels
end

const latex_table_column_labels = Dict{Any,String}(
        TranslationHasher => "No. of connected clusters",
        IsomorphicHasher => "No. of isomorphic clusters",
        SymmetricHasher => "No. of symmetric clusters",
)

_cluster_set_column_label(cluster_set::ClusterSet{C,H}) where {C,H} =
        get(latex_table_column_labels, H, "")

"""
    print_latex_table([io], expansion, cluster_sets, max_order)

Prints a LaTeX table akin to those found in scientific papers that utilize NLCE.
"""
function print_latex_table(io::IO, e::AbstractExpansion, cluster_sets::Vector{<:AbstractClusterSet}, max_order::Int)
        data, cluster_set_labels = _expansion_table_data(e, cluster_sets, max_order)
        cluster_set_headers = [LatexCell(l) for l in cluster_set_labels]
	column_labels = [vcat(["Order"], cluster_set_headers, [LatexCell("\$\\sum L(c)\$"), LatexCell("\$\\sum S(c)\$"), LatexCell("\$\\sum F(c)\$")])]
        style = LatexTableStyle(first_line_column_label=String[], column_label=String[])
        table_format = LatexTableFormat(; @latex__no_vertical_lines)
        pretty_table_latex_backend(io, data; column_labels=column_labels, style=style, table_format=table_format)
end

"""
    print_html_table([io], expansion, cluster_sets, max_order)

Prints an HTML table akin to those found in scientific papers that utilize NLCE.
"""
function print_html_table(io::IO, e::AbstractExpansion, cluster_sets::Vector{<:AbstractClusterSet}, max_order::Int)
        data, cluster_set_labels = _expansion_table_data(e, cluster_sets, max_order)
        column_labels = [vcat(["Order"], cluster_set_labels, ["sum L(c)", "sum |subgraphs|", "sum F(c)"])]
        pretty_table_html_backend(io, data; column_labels=column_labels)
end

"""
    print_ascii_table([io], expansion, cluster_sets, max_order)

Prints an ASCII table akin to those found in scientific papers that utilize NLCE.
"""
function print_ascii_table(io::IO, e::AbstractExpansion, cluster_sets::Vector{<:AbstractClusterSet}, max_order::Int)
        data, cluster_set_labels = _expansion_table_data(e, cluster_sets, max_order)
        column_labels = [vcat(["Order"], cluster_set_labels, ["sum L(c)", "sum |subgraphs|", "sum F(c)"])]
        pretty_table(io, data; column_labels=column_labels)
end

print_latex_table(e::AbstractExpansion, cluster_sets::Vector{<:AbstractClusterSet}, max_order::Int) =
        print_latex_table(stdout, e, cluster_sets, max_order)

print_html_table(e::AbstractExpansion, cluster_sets::Vector{<:AbstractClusterSet}, max_order::Int) =
        print_html_table(stdout, e, cluster_sets, max_order)

print_ascii_table(e::AbstractExpansion, cluster_sets::Vector{<:AbstractClusterSet}, max_order::Int) =
        print_ascii_table(stdout, e, cluster_sets, max_order)
