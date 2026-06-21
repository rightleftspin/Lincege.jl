"""
    ImportedCluster

A single cluster deserialized from a JSON file produced by `write_to_json`. To be used in subsequent physics algorithms
"""
struct ImportedCluster
        cluster_hash::UInt64
        order::Int
        n_sites::Int
        coordinates::Vector{Vector{Float64}}
        site_colors::Vector{Int}
        bonds::Vector{Tuple{Int,Int,Int}}
        weights::Vector{Float64}
end

"""
    import_from_json(filepath) -> Vector{ImportedCluster}

Read a JSON file produced by `write_to_json` and return a vector of
`ImportedCluster`, one per cluster entry in the file.
"""
function import_from_json(filepath::String)::Vector{ImportedCluster}
        data = JSON.parsefile(filepath)
        clusters = Vector{ImportedCluster}(undef, length(data))

        for (i, d) in enumerate(data)
                cluster_hash = parse(UInt64, d["cluster_hash"])
                order = d["order"]::Int
                n_sites = d["n_sites"]::Int
                coordinates = [Float64.(coord) for coord in d["coordinates"]]
                site_colors = Int.(d["site_colors"])
                bonds = Tuple{Int,Int,Int}[(b[1], b[2], b[3]) for b in d["bonds"]]
                weights = Float64.(d["weights"])

                clusters[i] = ImportedCluster(cluster_hash, order, n_sites, coordinates, site_colors, bonds, weights)
        end

        clusters
end

"""
    import_from_json_by_order(filepath) -> Vector{Vector{ImportedCluster}}

Like `import_from_json`, but returns clusters grouped by order. The returned
vector is indexed by order: `result[order]` holds all clusters of that order.
"""
function import_from_json_by_order(filepath::String)::Vector{Vector{ImportedCluster}}
        clusters = import_from_json(filepath)
        max_order = maximum(c.order for c in clusters)
        by_order = [Vector{ImportedCluster}() for _ in 1:max_order]
        for c in clusters
                push!(by_order[c.order], c)
        end
        by_order
end
