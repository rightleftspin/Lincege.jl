# Kagome Lattice

This example computes the linked cluster expansion on a shastry-sutherland up to
order 3.

## Setup

Define the unit cell of the Shastry-Sutherland lattice:

```@example shastry-sutherland
using Lincege

shsu_basis = [[0.0, 0.0], [(sqrt(3) + 1) / (2 * sqrt(2)), (sqrt(3) - 1) / (2 * sqrt(2))], [sqrt(2) / 2, sqrt(6) / 2], [(sqrt(3) + 3) / (2 * sqrt(2)), (sqrt(3) + 1) / (2 * sqrt(2))]]
shsu_pvecs = [[(sqrt(3) + 1) / sqrt(2), 0.0], [0.0, (sqrt(3) + 1) / sqrt(2)]]
shsu_bonds = [Bond(1, 2, [0, 0], 1), Bond(2, 3, [0, 0], 1), Bond(1, 4, [-1, 0], 1), Bond(3, 4, [-1, 0], 1),
        Bond(3, 4, [0, 0], 2), Bond(2, 3, [0, -1], 2), Bond(1, 2, [-1, 0], 2), Bond(1, 4, [-1, -1], 2),
        Bond(2, 4, [-1, -1], 3), Bond(1, 3, [-1, 0], 3)],
shsu_uc = UnitCell(shsu_basis, shsu_pvecs, shsu_bonds, [1, 2, 3, 4])
```
One of the important aspects of this lattice is coloring each basis site with a different color and the three distinct bond type which are represented with the last digit in the Bond() function. In this example of the shastry-sutherland lattice we have 4 different colors for the basis and three different types of bonds.

## Imaging lattice and unit cell
Using the defined unit cell it is possible to image both a singular unit cell and how the lattice images into a larger grid. Imaging the lattice is critical to ensure that all the information in the creation of the unit cell is correct. In order for the imaging functions to work you must first import the Plots library.

```@example shastry-sutherland
using Plots
Lincege.image_unit_cell(shsu_uc)
Lincege.image_lattice(shsu_uc)
```
## Building the lattice and clusters

```@example shastry-sutherland
m_order = 3
lattice = SiteExpansionLattice(m_order, shsu_uc)

trans_clusters = TranslationClusterSet(lattice)
clusters_from_lattice!(trans_clusters, lattice)

iso_clusters = IsomorphicClusterSet(lattice)
clusters_from_clusters!(iso_clusters, trans_clusters)
```

## Computing the expansion

```@example shastry-sutherland
expansion = Expansion(iso_clusters, lattice, m_order)
summation!(expansion, m_order)
expansion.weights
```

## Writing to JSON

```julia
write_to_json(expansion, lattice, iso_clusters, "shastry-sutherland_lattice.json")
```
