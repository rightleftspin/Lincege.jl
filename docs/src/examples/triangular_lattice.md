# Triangular Lattice

This example computes the linked cluster expansion on a triangular lattice up to
order 3.

## Setup

Define the triangular unit cell:

```@example kagome
using Lincege

triangular_basis = [[0.0, 0.0]]
triangular_pvecs = [[1.0, 0.0], [1.0/2, sqrt(3)/2]]
triangular_bonds = [Bond(1,1,[1,0],1),Bond(1,1,[0,1],1),Bond(1,1,[-1,1],1)]
triangular_uc = UnitCell(triangular_basis, triangular_pvecs, triangular_bonds, [1])
```

## Building the lattice and clusters

```@example kagome
m_order = 3
lattice = SiteExpansionLattice(m_order, triangular_uc)

trans_clusters = TranslationClusterSet(lattice)
clusters_from_lattice!(trans_clusters, lattice)

iso_clusters = IsomorphicClusterSet(lattice)
clusters_from_clusters!(iso_clusters, trans_clusters)
```

## Computing the expansion

```@example kagome
expansion = Expansion(iso_clusters, lattice, m_order)
summation!(expansion, m_order)
expansion.weights
```

## Writing to JSON

```julia
write_to_json(expansion, lattice, iso_clusters, "triangular_lattice.json")
```
