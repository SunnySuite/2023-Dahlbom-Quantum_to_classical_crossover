using DrWatson
@quickactivate "FiniteTemperatureFeI2"

using Sunny, GLMakie
include(srcdir("model.jl"))
include(srcdir("structure_factor_utils.jl"))


dims = (4, 4, 4)
sim_params = (;
    Δt_therm = 0.004,
    dur_therm = 12.5,
    dur_decorr = 5.0,
    Δt = 0.025,
    ωmax = 10.0,
    nω = 100,
    nsamples = 10,
    λ = 0.1,
)

sys, cryst = FeI2_sys_and_cryst((4, 4, 4))

κ = 1.070
kT = 0.1
@time estimate_sum(sys, κ, kT, sim_params)