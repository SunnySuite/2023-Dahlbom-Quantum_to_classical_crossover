using DrWatson
@quickactivate "FiniteTemperatureFeI2"

using Sunny
include(srcdir("model.jl"))
include(srcdir("sum_rule_utils.jl"))
include(srcdir("decorrelation_utils.jl"))

# After executing the first script, you should have some estimated κ values
# for a range of temperatures. In this script we will show how to estimate
# the dynamical spin structure factor at one of the temperatures for which
# a κ is available. We begin be loading the earlier results.
@unpack κs, kTs = load(datadir("kappas", "kappas.jld2"))

# In the second script, we found that T_N is approximately 3 K.
T_N = 3.05  # Neel temperature in K

# Choose one of the available temperatures above this value.
kT = kTs[10]
κ = κs[10]

# Our model is specified in terms of meV.
kT_meV = kT*Sunny.meV_per_K

# Set up the system parameters 
dims = (12, 12, 4)  # Manuscript using (24, 24, 8)
gs = 1              # Select one of the three ground states
seed = 101          # Seed for RNG

# Simulation parameters
Δt_therm = 0.004    # Step size for Langevin integrator
dur_therm = 10.0    # Safe thermalization time
λ = 0.1             # Phenomenological coupling to thermal bath
Δt = 0.025          # Integrator step size for dissipationless trajectories
nsamples = 10       # Number of dynamical trajectories to collect for estimating S(𝐪,ω).
                    # The manuscript used 1200.
ωmax = 10.0         # Maximum energy to resolve.
nω = 200            # Number of energy bins. Manuscript used 800.

# Estimate S(𝐪, ω), both with and without the renormalization. First set up
# sampled correlations object.
sys, cryst = FeI2_sys_and_cryst(dims; gs, seed)
sc = dynamical_correlations(sys; Δt, nω, ωmax)

# Thermalize the system.
thermalize!(sys, kT_meV, Δt_therm, dur_therm)

# Collect samples, applying the renormalization to the deterministic dynamics
# only.
@time for _ in 1:nsamples
    decorrelate!(sys, kT_meV, Δt_therm)  # Decorrelate the system using Langevin integration
    renormalize_system!(sys, κ)          # Turn on renormalization
    add_sample!(sc, sys)                 # Run a trajectory and calculate correlations
    renormalize_system!(sys, 1.0)        # Turn off renormalization so Langevin integration is unaffected
end

# Save the results.
data = DrWatson.@strdict sc sys cryst kT_meV
wsave(datadir("dssf", "renormalized.jld2"), data)

# And, for comparison, we'll calculate the the same thing without applying the
# renormalization. 
sys, cryst = FeI2_sys_and_cryst(dims; gs, seed)
sc = dynamical_correlations(sys; Δt, nω, ωmax)
thermalize!(sys, kT_meV, Δt_therm, dur_therm)
@time for _ in 1:nsamples
    decorrelate!(sys, kT_meV, Δt_therm) 
    add_sample!(sc, sys)
end

# Save the results.
data = DrWatson.@strdict sc sys cryst kT_meV
wsave(datadir("dssf", "unrenormalized.jld2"), data)