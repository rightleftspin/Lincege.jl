# Lincege.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://rightleftspin.github.io/Lincege.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://rightleftspin.github.io/Lincege.jl/dev/)
[![Build Status](https://github.com/rightleftspin/Lincege.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/rightleftspin/Lincege.jl/actions/workflows/CI.yml?query=branch%3Amain)

Lincege.jl is a julia package for computing LINked Cluster Expansions on a
GEneral geometry (LINCEGE). Utilizing Lincege.jl, you can easily generate
clusters and reduced lattice constants for any lattice in order to perform NLCE.

## Installation

In the Julia REPL type `]` to switch to package mode, then type the command:

```julia
pkg> add Lincege
```

## Try Lincege.jl

Generate square-lattice clusters containing up to three sites, grouping
them by graph isomorphism (topological clusters):

```julia
using Lincege

basis = [[0.0, 0.0]]
primitive_vectors = [[1.0, 0.0], [0.0, 1.0]]
bonds = [Bond(1, 1, [1, 0], 1), Bond(1, 1, [0, 1], 1)]
unit_cell = UnitCell(basis, primitive_vectors, bonds, [1])

max_order = 3
lattice = SiteExpansionLattice(max_order, unit_cell)

translation_clusters = TranslationClusterSet(lattice)
clusters_from_lattice!(translation_clusters, lattice)

isomorphic_clusters = IsomorphicClusterSet(lattice)
reduce_clusters!(isomorphic_clusters, translation_clusters)
```

To try Lincege.jl for an example use case, try the
[square lattice example](https://rightleftspin.github.io/Lincege.jl/dev/examples/square_lattice.html)
in the docs. For a more complicated site expansion example, try the
[kagome lattice example](https://rightleftspin.github.io/Lincege.jl/dev/examples/kagome_lattice.html).

Example scripts are in [`examples/`](examples/).
See the corresponding [README](examples/README.md) for run instructions.

## Visualizing a lattice

Load both `Plots` and `ColorSchemes` to activate the plotting extension. Run the 
following code to visualize the unit cell and corresponding lattice.

```julia
using Lincege, Plots, ColorSchemes

basis = [[0.0, 0.0]]
primitive_vectors = [[1.0, 0.0], [0.0, 1.0]]
bonds = [Bond(1, 1, [1, 0], 1), Bond(1, 1, [0, 1], 1)]
unit_cell = UnitCell(basis, primitive_vectors, bonds, [1])

image_unit_cell(unit_cell)
image_lattice(unit_cell)
```

## Contact and Citation

To contact the authors of this package, please email either Pranav Seetharaman
at [pjseetha@uwaterloo.ca](mailto:pjseetha@uwaterloo.ca) or Professor Ehsan
Khatami at [ehsan.khatami@sjsu.edu](mailto:ehsan.khatami@sjsu.edu).

If Lincege.jl helped you in your research, please cite us using the following
citation:

```bibtex
@misc{Lincege.jl,
	author  = {Seetharaman, Pranav and contributors},
	title   = {Lincege.jl},
	url     = {https://github.com/rightleftspin/Lincege.jl},
	version = {v1.0.0},
	year    = {2026},
	month   = {04}
}
```

Funded by the U.S. Department of Energy: Grant Number DE-SC0022311
