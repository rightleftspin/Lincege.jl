"""
    Cluster{
        tag_type<:Integer,
        mult_type<:Real
    }

A graph type with custom vertex labels containing vertex-, edge- and graph-level metadata.

# Fields
- `graph::MetaGraph`: underlying graph with vertex codes of type `Code`, includes details about the 
"""

struct Cluster{
    tag_type<:Integer,
    mult_type<:Real
}
    graph::MetaGraph
    subclusters::Vector{Cluster}
    tags::Vector{tag_type}
    multiplicity::mult_type
end

## Constructors, first the general constructor that should be used only if you know what you're doing
## or have some special case, then the more user friendly one based on the UnitCell, Bravais Lattice 
## and Neighbor Number

"""
    Cluster(cells, cell_colorings, edges, edge_info)

A constructor for the cluster struct for any lattice, finite or infinite, with specified
cells, which are vertex labels for the underlying graph. These cells have color, edges 
and corresponding edge weights which are in Vectors indexed in the same order as cells.

# Examples
```julia-repl
julia> Cluster([Cell([[0, 0]]), Cell([[1, 0]]), Cell([[0, 1]])], [1, 1, 1], [[2, 3], [1, 3], [1, 2]], [[[1], [1]], [[1], [1]], [[1], [1]]])
# TODO: LIST HOW THIS WORKS!
```
"""

function Cluster(
        cells::Vector{Cell},
        cell_colorings::Vector{Integer}, 
        edges::Vector{Vector{Integer}}, 
        edge_info::Vector{Vector{Vector{Integer}}},
    )
    cluster_graph = MetaGraph(
                              SimpleDiGraph(), 
                              label_type=Cell, 
                              vertex_data_type=Integer, 
                              edge_data_type=Vector{Integer}
                             )
    for i, cell in enumerate(cells)
        if !haskey(cluster_graph, cell)
            cluster_graph[cell] = cell_colorings[i]
        for j, neighbor in enumerate(edges[i])
            if haskey(cluster_graph, neighbor)
                cluster_graph[cell, cells[neighbor]] = edge_info[i][j]
            else
                cluster_graph[cells[neighbor]] = cell_colorings[neighbor]
                cluster_graph[cell, cells[neighbor]] = edge_info[i][j]
            end
        end
    end
    
    Cluster(
            cluster_graph,
            [],
            [],
            0
           )
end

"""
    Cluster(unit_cell, bravais_lattice, neighborhood, max_order, expand)

A constructor for the cluster struct for regular lattices that start with a specific 
bravais_lattice and unit_cell, tiles them up to max_order, make sure to specify a
max_order that is much higher than nessecary to ensure that every cluster is generated
in the infinite case, a large graph doesn't affect performance too much here.

The edges are connected by finding all points that are within the proper neighborhood 
and connecting them with edge weights populated both by direction and bond "length"

The unit cell may be expanded or not expanded depending on the expand flag, on default,
each unit cell is expanded till the individual cells only contain one point, upon which 
multiple center points are returned according to the number of sites in a unit cell.

Important Notes: 
    1) Here we use AlgebraicNumbers because exact numbers are nessecary for 
    the bravais lattice and neighborhood methods to work. We cannot accurately 
    determine which elements are within which neighborhood if we cannot perfectly 
    measure distance.

    2) Please make sure that the unit cell is within the bounds of the bravais lattice,
    or weird self intersection results might happen, ensure that you check your lattice
    before running intense computations. 

"""
function Cluster(
        unit_cell::Cell,
        bravais_lattice::Vector{Vector{AlgebraicNumber}}, 
        neighborhood::Int, 
        max_order::Int, 
        expand::Bool=true
    )

    dimensions = size(bravais_lattice)[1]
    coordinates::Vector{Vector{T}} = []
    # Cover each element in the unit cell
    for basis_elem in coordinates(unit_cell)
        # Generate all integers for coordinates in that basis
        for int_coord in 0:((max_order ^ dimensions) - 1)
            coordinate = []
            int_coord_temp = divrem(int_coord, max_order)
            # Find each coordinate 
            for dimension in 1:dimensions
                pushfirst!(coordinate, int_coord_temp[2])
                int_coord_temp = divrem(int_coord_temp[1], max_order)
            end
            # Shift according to the bravais lattice and basis
            push!(coordinates, sum(coordinate .* bravais_lattice) .+ basis_elem)
        end
    end

end

## Base extensions
function Base.show(io::IO, cluster::Cluster)
    graph_info = "Number of Vertices: $(nv(cluster))\n"
    graph_info *= "Number of Subgraphs: $(size(cluster.subclusters)[1])\n"
    graph_info *= "Multiplicity: $(cluster.multiplicity)\n"
    for label in labels(cluster)
        # Add the label information
        graph_info *= "Label: $(label)\n"
        # Add edges and their corresponding edge weights
        graph_info *= "\tEdges: $(size(outneighbor_labels(label)[1]))\n"
        for neighbor in outneighbor_labels(label)
            graph_info *= "\t\t$(neighbor) weights:$(cluster[label, neighbor])\n"
        end
    end
    # Print the total string with all the information so far.
    print(io, graph_info)
    return nothing
end

## Methods

#function nv()
#function allneighbors()
#function outneighbors()
#function inneighbors()
#function neighbor_labels()
#function outneighbor_labels()
#function inneighbor_labels()

## Link between graph codes and metagraph labels, rewritten here for the Cluster struct

"""
    code_for(cluster::Cluster, label)

Find the vertex code (or index) associated with label `label`.

This can be useful to pass to methods inherited from `Graphs`. Note, however, that vertex codes can be reassigned after vertex deletion.
"""
function code_for(cluster::Cluster, label)
    return cluster.graph.vertex_properties[label][1]
end

"""
    label_for(cluster::Cluster, code)

Find the label associated with code `code`.

This can be useful to interpret the results of methods inherited from `Graphs`. Note, however, that vertex codes can be reassigned after vertex deletion.
"""
function label_for(cluster::Cluster, code::Integer)
    return cluster.graph.vertex_labels[code]
end
