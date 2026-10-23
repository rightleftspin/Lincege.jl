# Finite Square Lattice with Periodic Boundaries

This example computes a cluster expansion on a finite 4×4 square lattice. Set `periodic=true`
to connect opposite edges in both directions. The lattice has 16 sites, each
with four neighbors. Here the expansion includes clusters of up to three sites.

## Setup

```@setup finite_square
using Lincege
if !@isdefined(output_dir)
    output_dir = joinpath(pkgdir(Lincege), "examples", "outputs")
end
mkpath(output_dir)
```

Define the one-site square unit cell:

```@example finite_square
using Lincege

square_basis = [[0.0, 0.0]]
square_pvecs = [[1.0, 0.0], [0.0, 1.0]]
square_bonds = [Bond(1, 1, [1, 0], 1), Bond(1, 1, [0, 1], 1)]
square_uc = UnitCell(square_basis, square_pvecs, square_bonds, [1])
```

## Building the lattice and clusters

Use `ConnectedClusterSet` to generate connected clusters on the finite lattice,
then reduce them by square symmetry or graph isomorphism:

```@example finite_square
finite_order = 3
finite_lattice = FiniteLattice((4, 4), square_uc, finite_order; periodic=true)

finite_connected = ConnectedClusterSet(finite_lattice)
clusters_from_lattice!(finite_connected, finite_lattice)

finite_symmetric = SymmetricClusterSet(finite_lattice, :Square)
reduce_clusters!(finite_symmetric, finite_connected)

finite_isomorphic = IsomorphicClusterSet(finite_lattice)
reduce_clusters!(finite_isomorphic, finite_connected)
```

## Computing the expansion and writing to JSON

Compute the expansion from the isomorphic clusters and export it for later.

```@example finite_square
finite_expansion = Expansion(finite_isomorphic, finite_lattice)
summation!(finite_expansion, finite_order)

finite_json_path = joinpath(output_dir, "square_lattice_periodic_finite_cluster_expansion.json")
write_to_json(finite_expansion, finite_lattice, finite_json_path)
```
