using DrWatson
@quickactivate "FiniteTemperatureFeI2"

using Sunny, GLMakie
include(srcdir("model.jl"))
include(srcdir("structure_factor_utils.jl"))


# General simulation parameters
dims = (4, 4, 4)  # Lattice size -- used (24, 24, 8) in paper
gs = 1            # Choose one of the three available ground states
seed = 1          # Seed for RNG
sys, cryst = FeI2_sys_and_cryst(dims; seed, gs)
sim_params = (;
    Δt_therm = 0.004,   # Langevin integrator step size, as in paper 
    λ = 0.1,            # Coupling to thermal bath, as in paper
    dur_therm = 12.5,   # Sufficient at most temperatures
    dur_decorr = 5.0,   # Sufficient at most temperatures
    Δt = 0.025,         # Sufficient for stability of dissipationless trajectories
    ωmax = 10.0,        # Maximum resolved energy
    nω = 100,           # Number of resolved energies -- 800 used in paper
    nsamples = 10,      # Number of sample trajectories
)


# Sum rule and κ parameters
ref = 16/3                     # Value of SU(3) quadratic casimir with chosen normalization convention
thresh = 0.05                  # Allowable deviation in estimated sum, relative to reference 
kTs = 10 .^ range(log10(0.1), log10(20.0), 10)  # 10 temperatures between 0.1 and 100.0, in meV
global_bounds = (1.0, 2.0)     # Smallest and largest κs (search space)


# Perform binary search to find κs at chosen temperatures
κs = zero(kTs)  
@time for (n, kT) in enumerate(kTs)
    println("kT = $kT")
    
    κ0 = n == 1 ? 1.0 : κs[n-1]  # Start search with κ=1 or previously found κ  
    bounds = (1.0, 2.0)          # Bounds on possible κ values
    total_weight = estimate_sum(sys, κ0, kT, sim_params; observables)

    @time while abs(total_weight - ref) > thresh
        total_weight = estimate_sum(sys, κ0, kT, sim_params; observables)
        println("\t κ=$κ0, sum=$total_weight")

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
    println("\t Final: κ=$κ0, sum=$total_weight")
    κs[n] = κ0
end

# For publication quality results, a larger system should be used and many more
# samples should be collected for each estimate of κ. Additionally, the same
# sampled initial conditions should be used for each kT to avoid the possibility
# of locking the binary search due to stochastic effects.

# Plot the results κs and save
scatter(kTs, κs; axis=(xscale=log10, xlabel="kT (meV)", ylabel="κ"))
data = DrWatson.@strdict kTs κs
wsave(datadir("kappas", "kappas.jld2"), data)