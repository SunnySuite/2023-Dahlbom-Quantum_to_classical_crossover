using DrWatson
@quickactivate "FiniteTemperatureFeI2"

using Sunny, GLMakie
include(srcdir("model.jl"))
include(srcdir("sum_rule_utils.jl"))

# # Finding κ(kT)

# The script demonstrates how to calculate moment renormalization factors, κ,
# described in the manuscript. These are constructed to ensure that a
# classically-calculated dynamical spin structure factor, S(𝐪,ω), satisfies the
# quantum sum rule. 
#
# ## κ search algorithm
#
# The approach is direct: for each T, select a κ and estimate ∫∫d𝐪dωS(𝐪,ω)
# using this κ. If the resulting sum is too large relative to the quantum sum
# rule, select a smaller κ; if the sum is too small, select a larger κ. To
# accelerate the search, we here write a function to perform a simple binary
# search.

function estimate_kappa(sys, kT, κ0, ref, global_bounds, sim_params; thresh=0.05, verbose=true)
    observables = observable_matrices() # Physical basis for SU(3)

    ## Initial estimate of spectral weight
    total_weight = estimate_sum(sys, κ0, kT, sim_params; observables)
    bounds = global_bounds

    ## Run binary search algorithm to find κ value that yields spectral weight
    ## sufficiently close to the reference
    @time while abs(total_weight - ref) > thresh
        total_weight = estimate_sum(sys, κ0, kT, sim_params; observables)
        println("a")
        verbose && println("\tκ0=$κ0, sum=$total_weight")
        if total_weight < ref
            if κ0 > global_bounds[2]
                κ0 = global_bounds[2]
                continue
            end
            bounds = (κ0, bounds[2])
            κ0 = (bounds[2] + κ0)/2
            continue
        elseif total_weight > ref
            if κ0 < global_bounds[1]
                κ0 = global_bounds[1]
                continue
            end
            bounds = (bounds[1], κ0)
            κ0 = (bounds[1] + κ0)/2
            continue
        end
    end

    return κ0
end


# ## Set up simulation parameters

dims = (4, 4, 4)  # Lattice size -- used (24, 24, 8) in paper
gs = 1            # Choose one of the three available ground states
seed = 1          # Seed for RNG
sys, cryst = FeI2_sys_and_cryst(dims; seed, gs)

sim_params = (;
    Δt_therm = 0.004,   # Langevin integrator step size (as in paper)
    λ = 0.1,            # Coupling to thermal bath (as in paper)
    dur_therm = 12.5,   # Sufficient for thermalization at most temperatures
    dur_decorr = 5.0,   # Sufficient for decorrelation at most temperatures
    Δt = 0.025,         # Sufficient for stability of dissipationless trajectories
    ωmax = 10.0,        # Maximum resolved energy
    nω = 100,           # Number of energy bins (800 used for manuscript)
    nsamples = 10,      # Number of sample trajectories (144 used for manuscript)
)


## Parameters for κ search 
global_bounds = (1.0, 2.0)  # Smallest and largest κs (search space)
thresh = 0.05               # Allowable deviation in estimated sum, relative to reference 
ref = 16/3                  # Quadratic Casimir of SU(3) with chosen normalization convention
kTs = 10 .^ range(log10(0.1), log10(30.0), 15)  # 15 temperatures between 0.1 and 10.0 in K

# ## Perform the search

## Estimate κs for chosen temperatures
κs = zero(kTs)  
for (n, kT) in enumerate(kTs)
    println("kT = $kT")
    κ0 = n == 1 ? 1.0 : κs[n-1]  # Start search with κ=1 or previously found κ  
    κs[n] = estimate_kappa(sys, kT, κ0, ref, global_bounds, sim_params; thresh)
end

# Plot the results,
scatter(kTs, κs; axis=(xscale=log10, xlabel="kT (meV)", ylabel="κ"))

# and save the data.

data = DrWatson.@strdict kTs κs
wsave(datadir("kappas", "kappas.jld2"), data)


# For publication quality results, a larger system should be used and many more
# samples should be collected for each estimate of κ(T). Additionally, the same
# sampled initial conditions should be used for each kT to avoid the possibility
# of locking the binary search due to stochastic effects -- e.g., the RNG should
# be reset to the same state before each trial in the binary search.