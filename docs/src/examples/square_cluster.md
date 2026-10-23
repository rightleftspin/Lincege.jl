# Square Lattice Square Expansion

This example computes the linked cluster expansion on the square lattice up to
order 3 using corner-sharing squares.

## Setup

```@setup square_cluster
using Lincege
if !@isdefined(output_dir)
    output_dir = joinpath(pkgdir(Lincege), "examples", "outputs")
end
mkpath(output_dir)
```

Define the unit cell with basis positions, primitive vectors, and
nearest-neighbour bonds:

```@example square_cluster
using Lincege

square_cluster_basis = [[[-1 / 2, -1 / 2], [-1 / 2, 1 / 2], [1 / 2, -1 / 2], [1 / 2, 1 / 2]]]
square_cluster_pvecs = [[1.0, 1.0], [1.0, -1.0]]
square_cluster_lbonds = [
        ExpansionBond([1, 1], [1, 2], [0, 0], 1),
        ExpansionBond([1, 1], [1, 3], [0, 0], 1),
        ExpansionBond([1, 2], [1, 4], [0, 0], 1),
        ExpansionBond([1, 3], [1, 4], [0, 0], 1),
]
square_cluster_ebonds = [
        Bond(1, 1, [1, 0], 1),
        Bond(1, 1, [0, 1], 1),
]
square_cluster_uc = ExpansionUnitCell(
        square_cluster_basis,
        square_cluster_pvecs,
        square_cluster_lbonds,
        square_cluster_ebonds,
        [[1, 1, 1, 1]]
)
```

Each expansion element is a four-site square. The `ExpansionBond` entries
connect physical sites, while the separate `Bond` entries connect neighboring
expansion elements. Equal site colors and bond types treat the sites and
interactions as equivalent.

## Imaging the unit cell and lattice

Load `Plots` and `ColorSchemes` to activate the plotting extension. The dashed
black lines outline the primitive cell; the solid lines are physical bonds.

```@example square_cluster
using Plots, ColorSchemes
image_unit_cell(square_cluster_uc)
```

The repeated lattice merges coincident shared sites. Enable the optional gray
dashed overlay to show connections between the centers of neighboring squares.
These connections describe the expansion geometry, rather than physical bonds.

```@example square_cluster
image_lattice(square_cluster_uc; show_expansion_bonds=true)
```

## Building the lattice and clusters

The following code generates clusters up to order 3, starting with
translational clusters and reducing to isomorphic clusters. Here order counts
expansion elements.

Use `WeakClusterExpansionLattice` because neighboring squares share corner
sites. 

```@example square_cluster
m_order = 3
lattice = WeakClusterExpansionLattice(m_order, square_cluster_uc)

trans_clusters = TranslationClusterSet(lattice)
clusters_from_lattice!(trans_clusters, lattice)

iso_clusters = IsomorphicClusterSet(lattice)
reduce_clusters!(iso_clusters, trans_clusters)
```

## Computing the expansion

The final NLCE expansion weights are computed recursively.

```@example square_cluster
expansion = Expansion(iso_clusters, lattice)
summation!(expansion, m_order)
```

## Writing to JSON

Lastly, write the results to a JSON file for future computations.

```julia
write_to_json(expansion, lattice, joinpath(output_dir, "square_lattice_weak_cluster_expansion.json"))
```
