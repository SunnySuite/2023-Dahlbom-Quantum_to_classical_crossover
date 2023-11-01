using DrWatson
@quickactivate "FiniteTemperatureFeI2"

using Sunny
include(srcdir("model.jl"))
include(srcdir("decorrelation_utils.jl"))

# After executing the first script, you should have some estimated κ values
# for a range of temperatures. In this script we will show how to estimate
# the dynamical spin structure factor at one of the temperatures for which
# a κ is available. We begin be loading the earlier results.

@unpack κs, kTs = load(datadir("kappas", "kappas.jld2"))

# In the second script, we found that T_N is approximately 3 K. We'll pick