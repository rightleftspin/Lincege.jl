```@meta
CurrentModule = Lincege
```

# Lincege

Documentation for [Lincege](https://github.com/rightleftspin/Lincege.jl), a
julia package for computing LINked Cluster Expansions on a GEneral geometry
(LINCEGE).

The general pipeline is:

1. Define a unit cell and lattice geometry.
2. Generate all unique clusters up to a desired order using a cluster set and hasher.
3. Build an `Expansion` from those clusters.
4. Call `summation!` to populate the NLCE weights.
5. Export results with `write_to_json`.

See the online documentation for worked examples on the square lattice, Kagome
lattice, and Pyrochlore unit cell.

```@autodocs
Modules = [Lincege]
```
