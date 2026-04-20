# One Fifth Depleted Lattice

This example computes the linked cluster expansion on a one fifth depleted lattice up to
order 3.

## Setup

Define the one-fifth-depleted unit cell:

[0,1,[0,0],0]
[0,2,[0,-1],0]
[1,3,[0,0],0]
[2,3,[0,1],0]
[1,2,[0,0],1]
[0,3,[-1,0],1]

```@example one-fifth-depleted
using Lincege
one_fifth_basis = [[0.0, 0.0], [0.0, 1.0], [0.0, 2.0],[1.0,1.0]]
one_fifth_pvecs = [[2.0, 1.0], [-1.0, 2.0]]
one_fifth_bonds = [Bond(1,2,[0,0],1),Bond(1,3,[0,-1],1),Bond(2,4,[0,0],1), Bond(3,4,[0,1],1),Bond(2,3,[0,0],2),Bond(1,4,[-1,0],2)]
one_fifth_uc = UnitCell(one_fifth_basis, one_fifth_pvecs, one_fifth_bonds, [1, 1, 1, 1])
```

## Building the lattice and clusters

```@example one-fifth-depleted
m_order = 3
lattice = SiteExpansionLattice(m_order, one_fifth_uc)

trans_clusters = TranslationClusterSet(lattice)
clusters_from_lattice!(trans_clusters, lattice)

iso_clusters = IsomorphicClusterSet(lattice)
clusters_from_clusters!(iso_clusters, trans_clusters)
```

## Computing the expansion

```@example one-fifth-depleted
expansion = Expansion(iso_clusters, lattice, m_order)
summation!(expansion, m_order)
expansion.weights
```

## Writing to JSON

```julia
write_to_json(expansion, lattice, iso_clusters, "one_fifth_lattice.json")
```
