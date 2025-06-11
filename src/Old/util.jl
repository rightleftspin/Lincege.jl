"""
General utility functions for a variety of parts of the NLCE process, these
will generally be math heavy functions that are used often
"""

"""
Generates the lattice by populating the basis around each site in the primitive lattice,
this adds it in terms of real space coordinates
"""
function add_basis_coords(basis, lattice)
        transpose(
                reduce(
                        hcat,
                        Iterators.flatten([
                                [site + elem for elem in basis] for site in eachcol(lattice)
                        ]),
                ),
        )
end

"""
Generates the lattice by populating the basis around each site in the primitive lattice,
this adds it in terms of sublattice coordinates
"""
function add_basis_sublattice(basis, unrotated_lattice)
        transpose(
                reduce(
                        hcat,
                        Iterators.flatten([
                                [[site..., i] for i = 1:length(basis)] for
                                site in eachcol(unrotated_lattice)
                        ]),
                ),
        )
end

"""
Generates every point from the primitive lattice vectors, in a cube of
side length = 2 * maximum_order + 1 centered at the origin
"""
function generate_primitive_lattice(primitive_vectors, max_order)

        unrotated_coords = generate_cartesian_coordinates(length(primitive_vectors), max_order)
        # rotate and stretch standard cartesian cube coordinates into primitive lattice
        (stack(primitive_vectors) * unrotated_coords, unrotated_coords)

end

"""
Generate every integer point in a cube of TWICE the given side length, centered around
the origin twice the side length is necessary to make sure a point exists at the origin
"""
function generate_cartesian_coordinates(dimension, half_side_length)
        # Forces the lattice to have a strict center point
        diameter = 2 * half_side_length + 1
        # Total number of coordinates for the entire lattice
        max_coords = diameter^dimension
        coords = repeat(transpose(0:(max_coords-1)), dimension)

        for dim = 0:(dimension-1)
                coords[dim+1, :] =
                        div.(coords[dim+1, :], diameter^dim) .% diameter .- half_side_length
        end

        coords
end

function adj_list_from_coords(coords, neighbor_distances)
        num_coords = size(coords, 1)
        adj_list::Vector{Vector{Int}} = []
        # Generates pairwise distances between all coordinates
        dists = pairwise(euclidean, coords, dims=1)

        for i = 1:num_coords
                temp_adj_list::Vector{Int} = []
                for j = 1:num_coords
                        for d in neighbor_distances
                                if (dists[i, j] ≈ d)
                                        append!(temp_adj_list, j)
                                end
                        end
                end
                push!(adj_list, temp_adj_list)
        end

        adj_list
end

function hashing_lattice_coords(
        real_space_coords,
        expansion_sublattice_coords,
        struct_per_basis,
        labels,
        translation_labels,
)
        sub_coords = []
        connection::Vector{Vector{Int}} = []
        rev_connection::Vector{Vector{Int}} = []
        all_labels = []
        all_translation_labels = []

        for (ind_exp, exp_coord, r_coord) in zip(
                1:length(eachrow(expansion_sublattice_coords)),
                eachrow(expansion_sublattice_coords),
                eachrow(real_space_coords),
        )
                temp_connection = []
                for (ind, sub_coord) in enumerate([
                        per_basis + r_coord for per_basis in struct_per_basis[exp_coord[end]]
                ])
                        sub_coord_loc = findfirst(≈(sub_coord), sub_coords)
                        if sub_coord_loc != nothing
                                append!(temp_connection, sub_coord_loc)
                                append!(rev_connection[sub_coord_loc], ind_exp)
                        else
                                push!(sub_coords, sub_coord)
                                push!(all_labels, labels[exp_coord[end]][ind])
                                push!(all_translation_labels, translation_labels[exp_coord[end]][ind])
                                append!(temp_connection, length(sub_coords))
                                push!(rev_connection, [ind_exp])
                        end
                end
                push!(connection, temp_connection)
        end

        (
                transpose(reduce(hcat, sub_coords)),
                connection,
                rev_connection,
                all_labels,
                all_translation_labels,
        )
end

"""
Converts adjacency list for a graph to an adjacency matrix
"""
function adj_list_to_adj_matrix(
        adj_list::AbstractVector{<:AbstractVector{<:Integer}},
        values::AbstractVector{<:AbstractVector{<:Integer}},
)

        # Find the number of vertices in the graph
        number_vertices = length(adj_list)
        # Initialize the empty adjacency matrix
        adj_matrix::Matrix{Int64} = zeros(number_vertices, number_vertices)

        # Loop over all the vertices in the adjacency list
        for (vertex, neighbors) in enumerate(adj_list)
                # Loop over all their neighbors
                for (neighbor_index, neighbor) in enumerate(neighbors)
                        # Add the corresponding value to the adjacency matrix
                        adj_matrix[vertex, neighbor] = values[vertex][neighbor_index]
                end
        end

        adj_matrix
end

function adj_list_to_adj_matrix(adj_list::AbstractVector{<:AbstractVector{<:Integer}})

        # Wrapper function for an adjacency list without weights. This function
        # will instead add a weight of 1 for every edge
        values = deepcopy(adj_list)
        for i = 1:length(values)
                values[i] .= 1
        end

        adj_list_to_adj_matrix(adj_list, values)
end

"""
Converts adjacency matrix for a graph to an adjacency list
"""
function adj_matrix_to_adj_list(adj_matrix::AbstractMatrix{<:Integer})

        # Find the number of vertices in the graph
        number_vertices = size(adj_matrix)[1]

        adj_list::Vector{Vector{Int64}} = []

        for i = 1:number_vertices
                neighbors = []
                for j = 1:number_vertices
                        if i != j
                                if adj_matrix[i, j] != 0
                                        append!(neighbors, j)
                                end
                        end
                end
                push!(adj_list, neighbors)
        end

        adj_list

end

"""
Converts adjacency matrix for a graph to an edge list"""
function adj_matrix_to_edge_list(adj_matrix::AbstractMatrix{<:Real})

        # Find the number of vertices in the graph
        number_vertices = size(adj_matrix)[1]

        edge_list::Vector{Vector{Int64}} = []

        for i = 1:number_vertices
                for j = i:number_vertices
                        if i != j
                                if adj_matrix[i, j] != 0
                                        push!(edge_list, [i, j, adj_matrix[i, j]])
                                end
                        end
                end
        end

        edge_list

end

"""
Below are more cluster related functions, they are not just purely mathematical functions.
"""
function translationally_invariant_clusters(
        lattice,
        start,
        max_order,
        single_site,
        per_site_factor,
)

        trans_invar_clusters = unique(
                c -> c[1],
                [
                        (Cluster(lattice, cluster_vertices), cluster_vertices) for
                        cluster_vertices in grow_lower(lattice, start, max_order)
                ],
        )

        if single_site
                append!(
                        trans_invar_clusters,
                        repeat(
                                [(
                                        Cluster(
                                                [Int64[]],
                                                Vector{Vector{Int}}(),
                                                [1; 1; 1;;;],
                                                vertex_labeled(trans_invar_clusters[1][1]),
                                                edge_labeled(trans_invar_clusters[1][1]),
                                        ),
                                        [1],
                                )],
                                per_site_factor,
                        ),
                )
        end

        (first.(trans_invar_clusters), last.(trans_invar_clusters))
end

function lattice_constants(hashing_fxn, per_site_factor::Integer, clusters, super_vertices)
        # (hash, (Cluster, Multiplicity, Permutation, super_vertices, subclusters(to be filled later)))
        cluster_info = Dict{UInt,Tuple{Cluster,Rational{Int},<:Any,<:Any,<:Any}}()
        add_mult_one =
                (cluster, mult, perm, svs, subs) ->
                        (cluster, mult + (1 // per_site_factor), perm, svs, subs)

        for (ind, (hash, permutation)) in enumerate(hashing_fxn.(clusters))

                cluster_info[hash] = add_mult_one(
                        get(
                                cluster_info,
                                hash,
                                (clusters[ind], 0, permutation, super_vertices[ind], []),
                        )...,
                )

        end

        cluster_info
end

function lattice_constants_only_info(
        hashing_fxn,
        per_site_factor::Integer,
        clusters,
        super_vertices,
)
        # (hash, (Cluster, Multiplicity, Permutation, super_vertices, subclusters(to be filled later)))
        cluster_info = Dict{UInt,Rational{Int}}()
        add_mult_one = (mult) -> (mult + (1 // per_site_factor))

        for (ind, (hash, permutation)) in enumerate(hashing_fxn.(clusters))

                cluster_info[hash] = add_mult_one(get(cluster_info, hash, 0)...)

        end

        cluster_info
end

function find_subclusters(cluster::Cluster, single_site::Bool)
        subclusters = []
        for order = 1:(nsv(cluster)-1)
                append!(
                        subclusters,
                        [
                                (Cluster(cluster, super_verts), super_verts) for
                                super_verts in grow_exact(cluster, super_vertices(cluster), order)
                        ],
                )
        end

        # TODO: Make this work for colored lattices as well
        if single_site
                append!(
                        subclusters,
                        repeat(
                                [(
                                        Cluster([Int64[]], Vector{Vector{Int}}(), [1; 1; 1;;;], false, false),
                                        [1],
                                )],
                                nv(cluster),
                        ),
                )
        end

        (first.(subclusters), last.(subclusters))
end

"""
Takes in a hashmap relating cluster hashes to their multiplicities, their
subclusters and their corresponding submultiplicities. Additionally, the
algorithm takes in the highest order to consider. Using this information,
the algorithm recursively generates a hashmap of the clusters and their
multiplicities. It will set the multiplicities of all clusters that are
higher than the given order to 0.
"""
function nlce_summation(clusters, order::Integer)
        cluster_weights = Dict()

        for (cluster_hash, (cluster, cluster_mult, _, _, _)) in clusters
                # Set clusters higher than the order to 0,
                if nsv(cluster) > order

                        cluster_weights[cluster_hash] = 0
                else
                        weights = _weight(clusters, cluster_hash)
                        map!(subcluster_weight -> cluster_mult * subcluster_weight, values(weights))
                        cluster_weights = mergewith(+, cluster_weights, weights)
                end
        end

        cluster_weights
end

"""
This is the recursive function that actually calculates the NLCE weights. It is
specified in many papers, I am using page 558 eqns 5 and 6 from the paper
"A short introduction to numerical linked-cluster expansions" by Tang, Khatami,
and Rigol. (https://arxiv.org/abs/1207.3366)

Inputs:
      clusters: Hashmap containing cluster hashes and their corresponding
      clusters, multiplicites, and subcluster information. Subcluster
      information is a hashmap of the subcluster hash and its associated
      multiplicity with reference to the cluster.

      cluster_hash: hash of the relevant cluster that we are finding
      the weights for

Output:
      Hashmap of cluster hashes and their corresponding multiplicity from
      the cluster specified in cluster_hash
"""
function _weight(clusters, cluster_hash::Integer)
        weight_dictionary = Dict([cluster_hash => 1 // 1])

        if nv(clusters[cluster_hash][1]) > 1
                for (subcluster_hash, subcluster_mult) in clusters[cluster_hash][5]
                        sub_weights = _weight(clusters, subcluster_hash)
                        map!(mult -> -1 * subcluster_mult * mult, values(sub_weights))
                        weight_dictionary = mergewith(+, weight_dictionary, sub_weights)
                end
        end

        weight_dictionary
end

"""
Resummation techniques and related functions. These resummation
techniques will take in a multidimensional array where the first
index is the property (energy, entropy, etc), the second index
is the independent variable (temperature, magnetic field, etc),
and the third index is over each order in the NLCE bare sum,
where each order corresponds to the partial sum up till that
order, from the weights returned by LINCEGE
"""
# TODO: Write the resummation functions more generally
function bincoeff(n, k)
        r = 1
        if (k > n)
                return 0
        end

        for d = 1:k
                r = r * n / d
                n -= 1
        end

        r
end

function break_properties(properties)
        broken_props = zero(properties)
        broken_props[:, :, 1] = properties[:, :, 1]

        for i = 2:size(properties, 3)
                broken_props[:, :, i] = properties[:, :, i] - properties[:, :, i-1]
        end

        broken_props
end

function euler_resummation(properties, start)
        broken_props = break_properties(properties[:, :, 1:(end-2)])
        partial_prop = broken_props[:, :, start:end]
        out = zero(partial_prop)
        for i = 1:size(partial_prop, 3)
                delta = zero(partial_prop[:, :, 1])
                for j = 1:i
                        coeff = bincoeff(i, j)
                        delta += (-1)^(j) * coeff .* abs.(partial_prop[:, :, j])
                end
                out[:, :, i] += (0.5^(i + 1)) .* delta
        end

        (sum(out, dims=3) + properties[:, :, start-1])
end

function eps_wynn(k, n, properties)
        if k == 0
                properties[:, :, n]
        elseif k == -1
                zero(properties[:, :, n])
        else
                first = eps_wynn(k - 2, n + 1, properties)
                second = eps_wynn(k - 1, n + 1, properties) .- eps_wynn(k - 1, n, properties)
                total = first + 1 ./ second
                total
        end
end

function wynn_resummation(properties, wynn_cycles)
        final_order = size(properties, 3)

        eps_wynn((2 * wynn_cycles), final_order - (2 * wynn_cycles), properties)
end

"""
These are writers that write a set of clusters to a file for future use. There are many styles
of writing to disk, but the standard will be JSON files
"""

function write_to_file(
        nlce_output::AbstractDict{<:Cluster,Vector{<:Real}},
        bundle::AbstractBundle,
        filename::AbstractString,
)
        clusters = []
        for (cluster, mults) in nlce_output
                cluster_dict = Dict()
                cluster_dict["NLCE Order"] = nsv(cluster)
                cluster_dict["Number of Sites"] = nv(cluster)
                cluster_dict["Site Labels"] = all_vertex_labels(cluster)
                # TODO: cluster_dict["coordinates"] = get_coordinates(bundle, cluster)
                cluster_dict["Weighted Bonds"] = weighted_edge_list(cluster)
                cluster_dict["Multiplicities"] = mults

                push!(clusters, cluster_dict)
        end

        open(filename, "w") do io
                JSON3.pretty(io, clusters)
        end
end

"""
Fortran format specifies for each cluster the number of bonds. Then it has each
bond listed below it in its own line, tab spaced. Between each cluster, there is
an empty line and at the end of the file is multiplicity for each cluster at
every order. For example, with a maximum order of 4 on a square lattice,
the file for the clusters of order 3 would look like this:

start of file >>>
2
1	2	1
1	3	1

0 0 6 -38
<<< end of file

Here there is only one cluster, it has 2 bonds. The bonds connect
vertices 1 and 2, and vertices 1 and 3 both with weight 1. The cluster
has overall multiplicity 0 for both order 1 and 2, multiplicity 6 for
order 3 and -38 for order 4.

A general example would look like this:

start of file >>>
number_of_bonds
vertex1 vertex2 bond_weight
vertex1 vertex3 bond_weight
...

number_of_bonds
vertex1 vertex2 bond_weight
vertex1 vertex3 bond_weight
...

multiplicity_1 multiplicity_2 ...
multiplicity_1 multiplicity_2 ...
...
<<< end of file
"""
function write_to_file_fortran(
        nlce_output::AbstractDict{<:Cluster,Vector{<:Real}},
        bundle::AbstractBundle,
        filename::AbstractString,
        max_order::Integer,
)

        nlce_files = [open(filename * "_$(i).txt", "w") for i = 1:max_order]
        sorted_clusters = sort(collect(keys(nlce_output)), by=nv)

        for cluster in sorted_clusters
                edges = weighted_edge_list(cluster)
                write(nlce_files[nv(cluster)], "$(length(edges))\n")
                for edge in edges
                        write(nlce_files[nv(cluster)], "$(join(edge, '\t'))\n")
                end
                write(nlce_files[nv(cluster)], "\n")
        end

        for cluster in sorted_clusters
                write(nlce_files[nv(cluster)], "$(join(nlce_output[cluster], ' '))\n")
        end

        close.(nlce_files)
end
