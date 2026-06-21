"""
    AbstractExpansionCluster

Abstract base type for an expansion cluster. Stores the subgraph information for a cluster, along with necessary information to write the cluster to disk.

Subtypes must implement:
- `subgraphs(cluster)` — subgraph hashes of each subgraph within the cluster
- `subtract_subcluster!(cluster, subcluster)` — subtract the weight (W_p(c)) of the given subcluster from the corresponding cluster 
"""
abstract type AbstractExpansionCluster end

lattice_constant(c::AbstractExpansionCluster) = _NI("lattice_constant")
subgraphs(c::AbstractExpansionCluster) = _NI("subgraphs")
subtract_subcluster!(c::AbstractExpansionCluster, sc::AbstractExpansionCluster) = _NI("subtract_subcluster!")

"""
    AbstractExpansion

Abstract base type for an NLCE expansion. Stores the subgraph relationships between
clusters and the weight matrix that `summation!` populates.

Subtypes must implement:
- `Base.getindex(e, cluster_hash)` — get the cluster corresponding to the `cluster_hash`
- `each_order(e, max_order)` — vector of clusters at each order up till `max_order`
"""
abstract type AbstractExpansion end

Base.getindex(e::AbstractExpansion, cluster_hash::UInt) = _NI("getindex")
each_order(e::AbstractExpansion, max_order::Int) = _NI("each_order")
order_offset(e::AbstractExpansion) = _NI("order_offset")
_expansion_table_data(e::AbstractExpansion, cs::Vector{<:AbstractClusterSet}, max_order::Int) = _NI("_expansion_table_data")

"""
    summation!(expansion, max_order)

Performs the recursive NLCE summation up till the max_order and populates the given expansion with the resultant weights
"""
summation!(e::AbstractExpansion, max_order::Int) = _summation!(e, max_order, order_offset(e))

function _summation!(e::AbstractExpansion, max_order::Int, order_offset::Int)
        for cluster_hashes in each_order(e, max_order + order_offset)
                for cluster_hash in cluster_hashes
                        cluster = e[cluster_hash]
                        for subcluster_hash in subgraphs(cluster)
                                subtract_subcluster!(cluster, e[subcluster_hash])
                        end
                end
        end
        e
end

include("ExpansionCluster.jl")
include("Expansion.jl")
include("util.jl")

const latex_table_column_labels = Dict{Any,String}(
        TranslationHasher => "No. of connected clusters",
        IsomorphicHasher => "No. of isomorphic clusters",
        SymmetricHasher => "No. of symmetric clusters",
)

_cs_column_label(cs::ClusterSet{C,H}) where {C,H} =
        get(latex_table_column_labels, H, "")

function print_latex_table(io::IO, e::AbstractExpansion, cs::Vector{<:AbstractClusterSet}, max_order::Int)
        data, cs_labels = _expansion_table_data(e, cs, max_order)
        cs_headers = [LatexCell(l) for l in cs_labels]
        column_labels = [vcat(["Order"], cs_headers, [LatexCell("\$\\sum L(c)\$"), LatexCell("\$\\sum |\\text{subgraphs}|\$"), LatexCell("\$\\sum \\F{c}\$")])]
        pretty_table_latex_backend(io, data; column_labels=column_labels)
end

function print_html_table(io::IO, e::AbstractExpansion, cs::Vector{<:AbstractClusterSet}, max_order::Int)
        data, cs_labels = _expansion_table_data(e, cs, max_order)
        column_labels = [vcat(["Order"], cs_labels, ["∑ L(c)", "∑ |subgraphs|", "∑ F(c)"])]
        pretty_table_html_backend(io, data; column_labels=column_labels)
end

function print_ascii_table(io::IO, e::AbstractExpansion, cs::Vector{<:AbstractClusterSet}, max_order::Int)
        data, cs_labels = _expansion_table_data(e, cs, max_order)
        column_labels = [vcat(["Order"], cs_labels, ["∑ L(c)", "∑ |subgraphs|", "∑ F(c)"])]
        pretty_table(io, data; column_labels=column_labels)
end

print_latex_table(e::AbstractExpansion, cs::Vector{<:AbstractClusterSet}, max_order::Int) =
        print_latex_table(stdout, e, cs, max_order)

print_html_table(e::AbstractExpansion, cs::Vector{<:AbstractClusterSet}, max_order::Int) =
        print_html_table(stdout, e, cs, max_order)

print_ascii_table(e::AbstractExpansion, cs::Vector{<:AbstractClusterSet}, max_order::Int) =
        print_ascii_table(stdout, e, cs, max_order)
