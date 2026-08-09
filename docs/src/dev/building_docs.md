# Building the Documentation for Lincege.jl

To build the docs locally, run `make.jl` from the `docs/` directory using the docs environment:

```bash
julia --project=. make.jl
```

Or equivalently from the repo root:

```bash
julia --project=docs docs/make.jl
```

Documenter.jl will write the rendered HTML to `docs/build/`. Open `docs/build/index.html` in a browser to view the result.

The warning:

```text
┌ Warning: Documenter could not auto-detect the building environment. Skipping deployment.
└ @ Documenter ~/.julia/packages/Documenter/AXNMp/src/deployconfig.jl:93
```

is expected — `deploydocs` only activates in CI (GitHub Actions) and is safe to ignore when building locally.
