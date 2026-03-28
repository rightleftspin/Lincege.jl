struct ClusterExpansion <: AbstractExpansion
        index_dictionary::Dict{UInt,Int}
        subgraphs::Vector{Vector{Int}}
        weights::Matrix{Float64}
        order_ids::Dict{Int,Vector{Int}}
end

function ClusterExpansion(clusters::AbstractClusterSet, lattice::AbstractClusterExpansionLattice, max_order::Int)
        index_dictionary = Dict{UInt,Int}()
        subgraphs = Vector{Vector{Int}}()
        order_ids = Dict{Int,Vector{Int}}(1 => [])
        weights = zeros(Float64, length(clusters), max_order + 1)

        for (i, cluster) in enumerate(sort(clusters))
                if haskey(order_ids, length(cluster) + 1)
                        push!(order_ids[length(cluster)+1], i)
                else
                        order_ids[length(cluster)+1] = [i]
                end
                index_dictionary[cluster.ghash] = i
                weights[i, length(cluster)+1] = cluster.lc
                temp_subgraphs = Int[]

                # Logic to add single site lvs
                for lv in connections(lattice)[cluster.evs]
                        gh = ghash(clusters, lv)

                        if !haskey(index_dictionary, gh)
                                ind = size(weights, 1) + 1
                                index_dictionary[gh] = ind
                                weights = vcat(weights, zeros(1, max_order + 1))
                                weights[ind, 1] = 1 / n_labels(lattice)
                                push!(order_ids[1], ind)
                        end

                        push!(temp_subgraphs, index_dictionary[gh])
                end

                for subgraph_evs in get_subgraphs(cluster, lattice)
                        push!(temp_subgraphs, index_dictionary[ghash(clusters, subgraph_evs)])
                end
                push!(subgraphs, temp_subgraphs)
        end

        ClusterExpansion(
                index_dictionary,
                subgraphs,
                weights,
                order_ids
        )
end

Base.getindex(e::ClusterExpansion, cluster_id::Int, order::Int) = getindex(e.weights, cluster_id, order)
Base.length(e::ClusterExpansion) = length(e.subgraphs)
order_ids(e::ClusterExpansion, order::Int) = e.order_ids[order]
get_subclusters(e::ClusterExpansion, cluster_id::Int) = e.subgraphs[cluster_id]
function add_array!(e::ClusterExpansion, order::Int, per_cluster::AbstractVector{Float64})
        @views e.weights[:, order] .+= per_cluster
end

summation!(e::ClusterExpansion, max_order::Int) = _summation!(e, max_order + 1)

function write_to_json(e::ClusterExpansion, lattice::AbstractClusterExpansionLattice, cs::AbstractClusterSet, filepath::String)
end
