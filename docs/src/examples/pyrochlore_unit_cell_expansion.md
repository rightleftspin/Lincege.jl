# Pyrochlore Lattice Unit Cell Expansion

This example computes the linked cluster expansion on the pyrochlore lattice up
to order 3 using four-site expansion cells.

## Setup

Define the unit cell with basis positions, primitive vectors, and
nearest-neighbour bonds:

```@example pyrochlore
using Lincege

pyro_exp_uc_basis = [
        [
                [1 / 2, 1 / 2, 1 / 2],
                [1 / 2, -1 / 2, -1 / 2],
                [-1 / 2, 1 / 2, -1 / 2],
                [-1 / 2, -1 / 2, 1 / 2],
        ],
]
pyro_exp_uc_pvecs = [[2.0, 2.0, 0.0], [2.0, 0.0, 2.0], [0.0, 2.0, 2.0]]
pyro_exp_uc_lbonds = [
        # Inside a Unit Cell
        ExpansionBond([1, 1], [1, 2], [0, 0, 0], 1),
        ExpansionBond([1, 1], [1, 3], [0, 0, 0], 1),
        ExpansionBond([1, 1], [1, 4], [0, 0, 0], 1),
        ExpansionBond([1, 2], [1, 3], [0, 0, 0], 1),
        ExpansionBond([1, 2], [1, 4], [0, 0, 0], 1),
        ExpansionBond([1, 3], [1, 4], [0, 0, 0], 1),
        # Between Unit Cells
        ExpansionBond([1, 1], [1, 2], [0, 0, 1], 1),
        ExpansionBond([1, 1], [1, 3], [0, 1, 0], 1),
        ExpansionBond([1, 1], [1, 4], [1, 0, 0], 1),
        ExpansionBond([1, 2], [1, 3], [0, 1, -1], 1),
        ExpansionBond([1, 2], [1, 4], [1, 0, -1], 1),
        ExpansionBond([1, 3], [1, 4], [1, -1, 0], 1),
]
pyro_exp_uc_ebonds = [
        Bond(1, 1, [1, 0, 0], 1),
        Bond(1, 1, [0, 1, 0], 1),
        Bond(1, 1, [0, 0, 1], 1),
        Bond(1, 1, [0, 1, -1], 1),
        Bond(1, 1, [1, 0, -1], 1),
        Bond(1, 1, [1, -1, 0], 1),
]
pyro_exp_uc_uc = ExpansionUnitCell(
        pyro_exp_uc_basis,
        pyro_exp_uc_pvecs,
        pyro_exp_uc_lbonds,
        pyro_exp_uc_ebonds,
        [[1, 1, 1, 1]]
)
```

Each expansion element contains four sites. The first six `ExpansionBond`
entries connect sites within the cell, the remaining six connect neighboring
cells. The separate `Bond` entries define connections between expansion
elements. All sites and physical bonds use label 1.

## Building the lattice and clusters

The following code generates clusters up to order 3, starting with
translational clusters and reducing to isomorphic clusters. Here order counts
expansion elements.

Use `StrongClusterExpansionLattice` because these expansion cells do not
share sites.

```@example pyrochlore
m_order = 3
lattice = StrongClusterExpansionLattice(m_order, pyro_exp_uc_uc)

trans_clusters = TranslationClusterSet(lattice)
clusters_from_lattice!(trans_clusters, lattice)

iso_clusters = IsomorphicClusterSet(lattice)
reduce_clusters!(iso_clusters, trans_clusters)
```

## Computing the expansion

The final NLCE expansion weights are computed recursively.

```@example pyrochlore
expansion = Expansion(iso_clusters, lattice)
summation!(expansion, m_order)
```

## Writing to JSON

Lastly, write the results to a JSON file for future computations.

```julia
write_to_json(expansion, lattice, "pyrochlore_lattice_uc_exp.json")
```
