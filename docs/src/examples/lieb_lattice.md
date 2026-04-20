# Lieb Lattice

This example computes the linked cluster expansion on a lieb lattice up to
order 3.

## Setup

Define the three-site lieb unit cell:

```@example lieb
using Lincege

lieb_basis = [[0.0, 0.0], [0.0, 1.0], [1, 0]]
lieb_pvecs = [[2.0, 0.0], [0.0, 2.0]]
lieb_bonds = [Bond(1,3,[1,0],1), Bond(1,2,[0,1],2)]
lieb_uc = UnitCell(lieb_basis, lieb_pvecs, lieb_bonds, [1, 1, 1])
```

## Building the lattice and clusters

```@example lieb
m_order = 3
lattice = SiteExpansionLattice(m_order, lieb_uc)

trans_clusters = TranslationClusterSet(lattice)
clusters_from_lattice!(trans_clusters, lattice)

iso_clusters = IsomorphicClusterSet(lattice)
clusters_from_clusters!(iso_clusters, trans_clusters)
```

## Computing the expansion

```@example lieb
expansion = Expansion(iso_clusters, lattice, m_order)
summation!(expansion, m_order)
expansion.weights
```

## Writing to JSON

```julia
write_to_json(expansion, lattice, iso_clusters, "lieb_lattice.json")
```
