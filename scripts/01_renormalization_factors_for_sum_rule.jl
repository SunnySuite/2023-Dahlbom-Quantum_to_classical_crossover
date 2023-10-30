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
κ = 1.0
kT = 0.1

@time estimate_sum(sys, κ, kT, sim_params; observables = observable_matrices())


kTs = range(0.1, 100.0, 3)
κs = zero(kTs)

ref = 16/3
thresh = 0.05
global_bounds = (1.0, 2.0)

for (n, kT) in enumerate(kTs)
    println("kT = $kT")
    
    κ0 = n == 1 ? 1.0 : κs[n-1]
    bounds = (1.0, 5.0)

    sum = estimate_sum(sys, κ0, kT, sim_params)

    @time while abs(sum - ref) > thresh
        sum = estimate_sum(sys, κ0, kT, sim_params)
        println("\t κ=$κ0, sum=$sum")
        if sum < ref
            if κ > global_bounds[2]
                κ0 = global_bounds[2]
                continue
            end
            bounds = (κ0, bounds[2])
            κ0 = (bounds[2] + κ0)/2
            continue
        elseif sum > ref
            if κ < global_bounds[1]
                κ0 = global_bounds[1]
                continue
            end
            bounds = (bounds[1], κ0)
            κ0 = (bounds[1] + κ0)/2
            continue
        end
    end

    κs[n] = κ0
end
