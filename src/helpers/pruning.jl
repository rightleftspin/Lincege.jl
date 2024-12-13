# Need to do
#
#       isomorphic_pruning:
#                       - account for edge weights

"""
These are pruning functions to be used inside step 2 of the algorithm,
these functions will always take in a cluster in some form and return
a hash of the cluster and a potential permutation of the cluster. If
the function returns nothing, this means that no permutation is necessary.
"""

"""
Takes a cluster and finds the translationally invariant
form of it.

Inputs: 
      cluster: takes a cluster struct, uses the direction
      weights matrix in the cluster

Output:
      Tuple of cluster hash and nothing
"""
function translational_pruning(cluster::AbstractNLCECluster)

    (
        hash(
            sum(
                weight -> 2^weight,
                direction_matrix(cluster)[sortperm(vertices(cluster)), :],
                dims = 2,
            ),
        ),
        nothing,
    )

end

"""
Finds the canonical ordering of the cluster using Nauty, this
function rearranges the cluster according to the permutation
used by Nauty

Inputs: 
      cluster: takes in a cluster struct, uses the
      edge_weighted_matrix inside it

Output:
      Tuple of cluster hash and permutation from nauty
"""
function isomorphic_pruning(cluster::AbstractNLCECluster)

    # Create an empty DenseNautyGraph
    nauty_graph =
        NautyGraph(edge_weighted_matrix(cluster), label(cluster, vertices(cluster)))

    # Canonize and find the corresponding permutation 
    permutation, _ = canonize!(nauty_graph)

    # Return the nauty hash and the permutation for the cluster
    # the slice is because the permutation will potentially
    # be longer than the initial graph, since edge weights
    # add extra vertices
    # Permutation goes backwards from the canonized graph to the
    # original graph
    return (ghash(nauty_graph), permutation[1:nv(cluster)])

end

"""
Takes a cluster and finds the translationally invariant
form of it. This does not return a hash.

Inputs: 
      cluster: takes a cluster struct, uses the direction
      weights matrix in the cluster

Output:
      Tuple of cluster hash and nothing
"""
function translational_pruning(cluster::AbstractNLCECluster)

            sum(
                weight -> 2^weight,
                direction_matrix(cluster)[sortperm(vertices(cluster)), :],
                dims = 2,
            )

end

function symmetric_pruning_pyrochlore(cluster::AbstractNLCECluster)
    L = nv(underlying_lattice(cluster))
    symmetric_variations = zeros(48, nv(cluster))
    for (i, n) in enumerate(vertices(cluster))
        
        symmetric_variations[:, i] = [
            L^2 * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) +
            L * (-L * fld(fld(n, L), L) + fld(n, L)) - L * fld(n, L) + n,
            -L^2 * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) + 2 * L^2 * fld(L, 2) - L * (-L * fld(fld(n, L), L) + fld(n, L)) +
            2 * L * fld(L, 2) - L * fld(n, L) + n,
            -L^2 * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) +
            2 * L^2 * fld(L, 2) +
            L * (-L * fld(fld(n, L), L) + fld(n, L)) +
            L * fld(n, L) - n + 2 * fld(L, 2),
            L^2 * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) -
            L * (-L * fld(fld(n, L), L) + fld(n, L)) +
            2 * L * fld(L, 2) +
            L * fld(n, L) - n + 2 * fld(L, 2),
            -L^2 * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) +
            2 * L^2 * fld(L, 2) +
            L * (-L * fld(n, L) + n) - L * fld(fld(n, L), L) + fld(n, L),
            L^2 * (-L * fld(n, L) + n) - L * (-L * fld(fld(n, L), L) + fld(n, L)) +
            2 * L * fld(L, 2) - L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L),
            L^2 * (-L * fld(fld(n, L), L) + fld(n, L)) +
            L * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) +
            L * fld(n, L) - n + 2 * fld(L, 2),
            -L^2 * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) + 2 * L^2 * fld(L, 2) - L * (-L * fld(n, L) + n) +
            2 * L * fld(L, 2) +
            L * fld(fld(n, L), L) +
            2 * fld(L, 2) - fld(n, L),
            -L^2 * (-L * fld(n, L) + n) + 2 * L^2 * fld(L, 2) -
            L * (-L * fld(fld(n, L), L) + fld(n, L)) +
            2 * L * fld(L, 2) +
            L * fld(fld(fld(n, L), L), L) +
            2 * fld(L, 2) - fld(fld(n, L), L),
            -L^2 * (-L * fld(fld(n, L), L) + fld(n, L)) + 2 * L^2 * fld(L, 2) -
            L * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) +
            2 * L * fld(L, 2) +
            L * fld(n, L) - n + 2 * fld(L, 2),
            L^2 * (-L * fld(fld(n, L), L) + fld(n, L)) + L * (-L * fld(n, L) + n) -
            L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L),
            -L^2 * (-L * fld(n, L) + n) + 2 * L^2 * fld(L, 2) -
            L * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) + 2 * L * fld(L, 2) -
            L * fld(fld(n, L), L) + fld(n, L),
            L^2 * (-L * fld(n, L) + n) -
            L * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) +
            2 * L * fld(L, 2) +
            L * fld(fld(n, L), L) +
            2 * fld(L, 2) - fld(n, L),
            -L^2 * (-L * fld(n, L) + n) +
            2 * L^2 * fld(L, 2) +
            L * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) +
            L * fld(fld(n, L), L) +
            2 * fld(L, 2) - fld(n, L),
            -L^2 * (-L * fld(fld(n, L), L) + fld(n, L)) +
            2 * L^2 * fld(L, 2) +
            L * (-L * fld(n, L) + n) +
            L * fld(fld(fld(n, L), L), L) +
            2 * fld(L, 2) - fld(fld(n, L), L),
            -L^2 * (-L * fld(fld(n, L), L) + fld(n, L)) + 2 * L^2 * fld(L, 2) -
            L * (-L * fld(n, L) + n) + 2 * L * fld(L, 2) - L * fld(fld(fld(n, L), L), L) +
            fld(fld(n, L), L),
            L^2 * (-L * fld(fld(n, L), L) + fld(n, L)) - L * (-L * fld(n, L) + n) +
            2 * L * fld(L, 2) +
            L * fld(fld(fld(n, L), L), L) +
            2 * fld(L, 2) - fld(fld(n, L), L),
            L^2 * (-L * fld(n, L) + n) +
            L * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) -
            L * fld(fld(n, L), L) + fld(n, L),
            L^2 * (-L * fld(fld(n, L), L) + fld(n, L)) -
            L * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) + 2 * L * fld(L, 2) -
            L * fld(n, L) + n,
            -L^2 * (-L * fld(n, L) + n) +
            2 * L^2 * fld(L, 2) +
            L * (-L * fld(fld(n, L), L) + fld(n, L)) - L * fld(fld(fld(n, L), L), L) +
            fld(fld(n, L), L),
            L^2 * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) +
            L * (-L * fld(n, L) + n) +
            L * fld(fld(n, L), L) +
            2 * fld(L, 2) - fld(n, L),
            -L^2 * (-L * fld(fld(n, L), L) + fld(n, L)) +
            2 * L^2 * fld(L, 2) +
            L * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) - L * fld(n, L) + n,
            L^2 * (-L * fld(n, L) + n) +
            L * (-L * fld(fld(n, L), L) + fld(n, L)) +
            L * fld(fld(fld(n, L), L), L) +
            2 * fld(L, 2) - fld(fld(n, L), L),
            L^2 * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) -
            L * (-L * fld(n, L) + n) + 2 * L * fld(L, 2) - L * fld(fld(n, L), L) +
            fld(n, L),
            -L^2 * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) + 2 * L^2 * fld(L, 2) - L * (-L * fld(fld(n, L), L) + fld(n, L)) +
            2 * L * fld(L, 2) +
            L * fld(n, L) - n + 2 * fld(L, 2),
            -L^2 * (-L * fld(fld(n, L), L) + fld(n, L)) +
            2 * L^2 * fld(L, 2) +
            L * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) +
            L * fld(n, L) - n + 2 * fld(L, 2),
            L^2 * (-L * fld(n, L) + n) - L * (-L * fld(fld(n, L), L) + fld(n, L)) +
            2 * L * fld(L, 2) +
            L * fld(fld(fld(n, L), L), L) +
            2 * fld(L, 2) - fld(fld(n, L), L),
            -L^2 * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) + 2 * L^2 * fld(L, 2) - L * (-L * fld(n, L) + n) + 2 * L * fld(L, 2) -
            L * fld(fld(n, L), L) + fld(n, L),
            L^2 * (-L * fld(fld(n, L), L) + fld(n, L)) -
            L * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) +
            2 * L * fld(L, 2) +
            L * fld(n, L) - n + 2 * fld(L, 2),
            -L^2 * (-L * fld(n, L) + n) + 2 * L^2 * fld(L, 2) -
            L * (-L * fld(fld(n, L), L) + fld(n, L)) + 2 * L * fld(L, 2) -
            L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L),
            -L^2 * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) +
            2 * L^2 * fld(L, 2) +
            L * (-L * fld(n, L) + n) +
            L * fld(fld(n, L), L) +
            2 * fld(L, 2) - fld(n, L),
            -L^2 * (-L * fld(fld(n, L), L) + fld(n, L)) + 2 * L^2 * fld(L, 2) -
            L * (-L * fld(n, L) + n) +
            2 * L * fld(L, 2) +
            L * fld(fld(fld(n, L), L), L) +
            2 * fld(L, 2) - fld(fld(n, L), L),
            L^2 * (-L * fld(n, L) + n) +
            L * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) +
            L * fld(fld(n, L), L) +
            2 * fld(L, 2) - fld(n, L),
            -L^2 * (-L * fld(n, L) + n) +
            2 * L^2 * fld(L, 2) +
            L * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) -
            L * fld(fld(n, L), L) + fld(n, L),
            L^2 * (-L * fld(n, L) + n) -
            L * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) + 2 * L * fld(L, 2) -
            L * fld(fld(n, L), L) + fld(n, L),
            L^2 * (-L * fld(fld(n, L), L) + fld(n, L)) - L * (-L * fld(n, L) + n) +
            2 * L * fld(L, 2) - L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L),
            L^2 * (-L * fld(fld(n, L), L) + fld(n, L)) +
            L * (-L * fld(n, L) + n) +
            L * fld(fld(fld(n, L), L), L) +
            2 * fld(L, 2) - fld(fld(n, L), L),
            -L^2 * (-L * fld(fld(n, L), L) + fld(n, L)) +
            2 * L^2 * fld(L, 2) +
            L * (-L * fld(n, L) + n) - L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L),
            -L^2 * (-L * fld(n, L) + n) + 2 * L^2 * fld(L, 2) -
            L * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) +
            2 * L * fld(L, 2) +
            L * fld(fld(n, L), L) +
            2 * fld(L, 2) - fld(n, L),
            L^2 * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) +
            L * (-L * fld(fld(n, L), L) + fld(n, L)) +
            L * fld(n, L) - n + 2 * fld(L, 2),
            L^2 * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) -
            L * (-L * fld(fld(n, L), L) + fld(n, L)) + 2 * L * fld(L, 2) - L * fld(n, L) +
            n,
            -L^2 * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) +
            2 * L^2 * fld(L, 2) +
            L * (-L * fld(fld(n, L), L) + fld(n, L)) - L * fld(n, L) + n,
            L^2 * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) -
            L * (-L * fld(n, L) + n) +
            2 * L * fld(L, 2) +
            L * fld(fld(n, L), L) +
            2 * fld(L, 2) - fld(n, L),
            -L^2 * (-L * fld(n, L) + n) +
            2 * L^2 * fld(L, 2) +
            L * (-L * fld(fld(n, L), L) + fld(n, L)) +
            L * fld(fld(fld(n, L), L), L) +
            2 * fld(L, 2) - fld(fld(n, L), L),
            -L^2 * (-L * fld(fld(n, L), L) + fld(n, L)) + 2 * L^2 * fld(L, 2) -
            L * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) + 2 * L * fld(L, 2) -
            L * fld(n, L) + n,
            L^2 * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) +
            L * (-L * fld(n, L) + n) - L * fld(fld(n, L), L) + fld(n, L),
            L^2 * (-L * fld(n, L) + n) + L * (-L * fld(fld(n, L), L) + fld(n, L)) -
            L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L),
            L^2 * (-L * fld(fld(n, L), L) + fld(n, L)) +
            L * (-L * fld(fld(fld(n, L), L), L) + fld(fld(n, L), L)) - L * fld(n, L) + n,
        ]
    end

    cluster_orbit = unique(sort(symmetric_variations, dims = 2), dims = 1)
    cluster_orbit_centered = cluster_orbit .- cluster_orbit[:, 1]

    (hash(sortslices(cluster_orbit_centered, dims=1)), nothing)
end

function symmetric_pruning_square(cluster::AbstractNLCECluster)

end

function symmetric_pruning_triangular(cluster::AbstractNLCECluster)

end
