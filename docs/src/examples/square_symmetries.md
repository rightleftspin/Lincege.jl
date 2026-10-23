# Square Lattice Symmetric Clusters

This example generates square-lattice clusters up to order 4 and compares
translation, point-group symmetry, and graph-isomorphism reductions.

## Setup

```@setup square_symmetries
using Lincege
if !@isdefined(output_dir)
    output_dir = joinpath(pkgdir(Lincege), "examples", "outputs")
end
mkpath(output_dir)
```

Define the square unit cell and bonds:

```@example square_symmetries
using Lincege

basis = [[0.0, 0.0]]
primitive_vectors = [[1.0, 0.0], [0.0, 1.0]]
bonds = [Bond(1, 1, [1, 0], 1), Bond(1, 1, [0, 1], 1)]
unit_cell = UnitCell(basis, primitive_vectors, bonds, [1])

max_order = 4
lattice = SiteExpansionLattice(max_order, unit_cell)
```

## Explicit symmetry matrices

These are the eight square lattice matrices which act on Cartesian position
vectors and form the square's point group.

Available geometry presets are `:Square`, `:Shastry_Sutherland`, `:Pyrochlore`,
`:Triangular`, and `:Kagome`; use `SymmetricClusterSet(lattice, :Square)` to
select a preset instead of supplying the matrices explicitly.

```@example square_symmetries
symmetries = Vector{Matrix{Float64}}([
        [1 0; 0 1],
        [0 1; -1 0],
        [-1 0; 0 -1],
        [0 -1; 1 0],
        [0 1; 1 0],
        [0 -1; -1 0],
        [1 0; 0 -1],
        [-1 0; 0 1],
])
```

## Generating and reducing clusters

First generate the translational clusters then reduce them into the
symmetric clusters and finally isomorphic clusters.

```@example square_symmetries
translation_clusters = TranslationClusterSet(lattice)
clusters_from_lattice!(translation_clusters, lattice)

symmetric_clusters = SymmetricClusterSet(lattice, symmetries)
reduce_clusters!(symmetric_clusters, translation_clusters)

isomorphic_clusters = IsomorphicClusterSet(lattice)
reduce_clusters!(isomorphic_clusters, symmetric_clusters)
```

## Computing the expansions

Build an `Expansion` for each reduction and compute its weights:

```@example square_symmetries
symmetric_expansion = Expansion(symmetric_clusters, lattice)
summation!(symmetric_expansion, max_order)

isomorphic_expansion = Expansion(isomorphic_clusters, lattice)
summation!(isomorphic_expansion, max_order)
```

## Relating symmetric and isomorphic clusters

Each isomorphic entry includes its weights and a `related_symmetric` list 
containing the symmetric cluster's weights and coordinates reordered to
match the canonical graph. 

```@example square_symmetries
json_path = joinpath(output_dir, "square_lattice_site_expansion_isomorphic_symmetric.json")
write_to_json(
        isomorphic_expansion, isomorphic_clusters,
        symmetric_expansion, symmetric_clusters,
        lattice, json_path,
)
```

The text export preserves every field of the joint JSON export:

```@example square_symmetries
txt_path = joinpath(output_dir, "square_lattice_site_expansion_isomorphic_symmetric.txt")
write_to_txt_file(
        isomorphic_expansion, isomorphic_clusters,
        symmetric_expansion, symmetric_clusters,
        lattice, txt_path,
)
```

The colon-separated records have two layouts:

```text
isomorphic : n_sites : cluster_hash : canonical_bonds : canonical_site_colors : weights : order
symmetric : isomorphic_hash : cluster_hash : permuted_coordinates : weights
```

Each isomorphic record is followed by its related symmetric records. The parent
hash links the two reductions; both retain their own hashes and weights.
Bonds are space-separated triples, coordinates use `(x,y)` tuples, and empty
fields represent empty vectors. Header lines beginning with `#` name the fields.
