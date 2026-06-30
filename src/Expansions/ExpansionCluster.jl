"""
    ExpansionCluster(cluster, lattice)

Cluster from a set of clusters that contains information necessary for summation! and writing to disk.
"""
struct ExpansionCluster <: AbstractExpansionCluster
        vertices::LatticeVertices{Int}
        lattice_constant::Float64
        subgraphs::Vector{UInt}
        weights::Dict{UInt,Float64}
end

function ExpansionCluster(cluster::AbstractCluster, clusters::AbstractClusterSet, lattice::SiteExpansionLattice)
        cluster_lattice_constant = lattice_constant(cluster)
        subgraphs = Vector{UInt}()

        for subgraph in get_subgraphs(cluster, lattice)
                push!(subgraphs, ghash(clusters, subgraph))
        end

        ExpansionCluster(cluster.vertices, cluster_lattice_constant, subgraphs, Dict{UInt,Float64}(cluster.ghash => 1.0))
end

function ExpansionCluster(cluster::AbstractCluster, clusters::AbstractClusterSet, lattice::StrongClusterExpansionLattice)
        cluster_lattice_constant = lattice_constant(cluster)
        lattice_vertices = connections(lattice)[cluster.vertices]
        subgraphs = Vector{UInt}()
        sizehint!(subgraphs, length(lattice_vertices) + length(cluster) - 1)

        for lv in lattice_vertices
                push!(subgraphs, ghash(clusters, LatticeVertices(lv)))
        end
        for subgraph in get_subgraphs(cluster, lattice)
                push!(subgraphs, ghash(clusters, subgraph))
        end

        ExpansionCluster(lattice_vertices, cluster_lattice_constant, subgraphs, Dict{UInt,Float64}(cluster.ghash => 1.0))
end

function ExpansionCluster(cluster::AbstractCluster, clusters::AbstractClusterSet, lattice::WeakClusterExpansionLattice)
        cluster_lattice_constant = lattice_constant(cluster)
        lattice_vertices = just_lattice_vertices(connections(lattice), cluster.vertices)
        subgraphs = Vector{UInt}()
        sizehint!(subgraphs, length(lattice_vertices) + length(cluster) - 1)

        for lv in lattice_vertices
                push!(subgraphs, ghash(clusters, LatticeVertices(lv)))
        end
        for subgraph in get_subgraphs(cluster, lattice)
                push!(subgraphs, ghash(clusters, subgraph))
        end

        ExpansionCluster(lattice_vertices, cluster_lattice_constant, subgraphs, Dict{UInt,Float64}(cluster.ghash => 1.0))
end

function ExpansionCluster(lv::Int, single_site_hash::UInt, n_single_site_clusters::Int)
        cluster_lattice_constant = 1 / n_single_site_clusters

        ExpansionCluster(LatticeVertices(lv), cluster_lattice_constant, UInt[], Dict{UInt,Float64}(single_site_hash => 1.0))
end

lattice_constant(cluster::ExpansionCluster) = cluster.lattice_constant
subgraphs(cluster::ExpansionCluster) = cluster.subgraphs
function subtract_subcluster!(cluster::ExpansionCluster, subcluster::ExpansionCluster)
        for (k, v) in subcluster.weights
                cluster.weights[k] = get(cluster.weights, k, 0.0) - v
        end
end
