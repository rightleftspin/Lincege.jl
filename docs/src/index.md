```@meta
CurrentModule = Lincege
```

# Lincege

This is the documentation for [Lincege](https://github.com/rightleftspin/Lincege.jl), a
julia package for computing LINked Cluster Expansions on a GEneral geometry (LINCEGE).

The general pipeline that LINCEGE follows is:
1. Define a unit cell and lattice geometry.
2. Generate all unique clusters up to a desired order using a cluster set and hasher.
3. Build an `Expansion` from those clusters.
4. Call `summation!` to populate the NLCE weights.
5. Export results with `write_to_json`.

Worked examples are shown in the side bar. A full computation from lattice information 
to computed observables for the square lattice Ising model is shown in 
[Square Lattice Ising Model](examples/square_lattice.md), which serves as a good starting point for new users.

For advanced users, see [API Reference](api.md) for the complete documentation.
