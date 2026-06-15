# Square Lattice

This example computes the linked cluster expansion on a square lattice up to
order 3.

## Setup

Define the unit cell with basis positions, primitive vectors, and
nearest-neighbour bonds:

```@example square
using Lincege

square_basis = [[0.0, 0.0]]
square_pvecs = [[1.0, 0.0], [0.0, 1.0]]
square_bonds = [Bond(1, 1, [1, 0], 1), Bond(1, 1, [0, 1], 1)]
square_uc = UnitCell(square_basis, square_pvecs, square_bonds, [1])
```

## Building the lattice and clusters

```@example square
m_order = 3
lattice = SiteExpansionLattice(m_order, square_uc)

trans_clusters = TranslationClusterSet(lattice)
clusters_from_lattice!(trans_clusters, lattice)

iso_clusters = IsomorphicClusterSet(lattice)
clusters_from_clusters!(iso_clusters, trans_clusters)

sym_clusters = SymmetricClusterSet(lattice, :Square)
clusters_from_clusters!(sym_clusters, trans_clusters)
```

## Computing the expansion

```@example square
expansion = Expansion(iso_clusters, lattice)
summation!(expansion, m_order)
```

## Printing Tables

There are three ways to print the standard tables one sees in papers regarding
NLCE.

`print_latex_table` outputs LaTeX source ready to paste into a paper:

```@example square
print_latex_table(expansion, [trans_clusters, iso_clusters, sym_clusters], m_order)
```

`print_html_table` renders the table in the browser:

```@example square
using Markdown
io = IOBuffer()
Lincege.print_html_table(io, expansion, [trans_clusters, iso_clusters, sym_clusters], m_order)
Markdown.HTML(String(take!(io)))
```

`print_ascii_table` prints a plain-text table:

```@example square
print_ascii_table(expansion, [trans_clusters, iso_clusters, sym_clusters], m_order)
```

## Writing to JSON

```julia
write_to_json(expansion, lattice, iso_clusters, "square_lattice.json")
```
