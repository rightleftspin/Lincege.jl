# Examples

Run all commands below from the repository root, the parent directory of `examples/`.

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
julia -t 4 --project=docs examples/square_lattice.jl
```
File outputs are saved in `examples/outputs/`, which is created automatically.
Note that the above command is run
with 4 threads (`-t 4`), for orders below 12 on most lattices, this does not
provide a large performance gain.
