# Examples

The `.jl` scripts in this directory are automatically extracted from the project
documentation using the following command.
```bash
julia --project=docs docs/example_scripts.jl
```

To run an example script, first instantiate the docs dependencies with
```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
```

Then run an example with:
```bash
julia --project=docs examples/square_lattice.jl
```
