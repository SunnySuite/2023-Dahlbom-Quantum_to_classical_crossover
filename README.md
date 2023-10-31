# Quantum-to-classical crossover in generalized spin systems -- the temperature-dependent spin dynamics of FeI2

The repository demonstrates the basic procedures for reproducing the calculations in the 2023 manuscript, ["Quantum-to-classical crossover in generalized spin systems -- the temperature-dependent spin dynamics of FeI2".](https://put.link.here.when.available)

## Notes on implementation and execution

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> FiniteTemperatureFeI2

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "2023-Dahlbom-Quantum_to_classical_crossover"
```
which auto-activate the project and enable local path handling from DrWatson.
