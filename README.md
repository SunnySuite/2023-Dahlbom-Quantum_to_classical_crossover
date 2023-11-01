# Quantum-to-classical crossover in generalized spin systems -- the temperature-dependent spin dynamics of FeI2

The repository contains the code necessary to implement the calculations in the 2023 manuscript, [Quantum-to-classical crossover in generalized spin systems -- the temperature-dependent spin dynamics of FeI2](https://arxiv.org/abs/2310.19905). It contains four Julia scripts using the [Sunny.jl](https://github.com/SunnySuite/Sunny.jl) library. These demonstrate:

1) How to calculate the renormalization factors described in the manuscript above
2) How to calculate Bragg intensities for the determination of the classical Néel temperature
3) The estimation of the dynamical spin structure factor using the renormalization factors
4) The extraction and presentation of intensities

The scripts present the minimal necessary machinery to perform a complete calculation. Parameters are selected so that execution time is minimal. Full reproduction of the figures in the paper can be achieved in most cases by increasing systems sizes, increasing energy resolution, examining a wider range of temperatures, and generating many more samples.

A short-form tutorial demonstrating the approach in even more abbreviated form is forthcoming.

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
