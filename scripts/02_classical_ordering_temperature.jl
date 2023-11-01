using DrWatson
@quickactivate "FiniteTemperatureFeI2"

using Sunny, Statistics, GLMakie
include(srcdir("model.jl"))
include(srcdir("decorrelation_utils.jl"))

# Simulation parameters
dur_therm = 15.0    # Conservative thermalization time for all temperatures
Δt = 0.004          # Step size sufficient to ensure stability
λ = 0.1             # Phenomenological coupling to thermal bath


# Set up a System
dims = (12, 12, 4)
gs = 1 seed = 101


# Estimate the Bragg intensities at a range of temperatures
kTs = 10 .^ range(log10(0.1), log10(10), 20) # Temperature range in Kelvin. 
                                             # Manuscript used 88 temperatures between 0.1 and 77 Kelvin.
nsamples = 10   # Manuscript used 3000 samples per temperature

μs = Float64[]
σs = Float64[]
for kT in kTs
    # Initialize a system in the ground state and prepare instant correlation calculation
    sys, cryst = FeI2_sys_and_cryst(dims; gs, seed)
    ic = instant_correlations(sys)
    formula = intensity_formula(ic, :trace)
    q_ordering = [q_gs[gs]]  # Ordering wave vector for chosen ground state

    # Thermalize the system at kT
    kT_meV = kT * Sunny.meV_per_K
    thermalize!(sys, kT_meV, Δt, dur_therm; λ)

    # Collect samples 
    bragg_samples = zeros(nsamples)
    println("Sample at kT=$kT")
    @time for n in 1:nsamples
        # Decorrelate the system and collect a sample
        decorrelate!(sys, kT_meV, Δt)
        add_sample!(ic, sys)

        # Evaluate the Bragg intensity, I(q_ord)
        bragg_samples[n] = instant_intensities_interpolated(ic, q_ordering, formula) |> only

        # Clear out the data in the correlation sampler so no cumulative average
        # is taken. We instead wish to collect the sample Bragg intensities one
        # at a time (as in the line above) so we can estimate the standard
        # error.
        ic.data .= 0.0
        ic.nsamples[1] = 0
    end

    # Save the mean and std of the Bragg intensities
    push!(μs, mean(bragg_samples))
    push!(σs, std(bragg_samples))
end


# Calculate numerical derivative of the Bragg intensities as a function of temperature.
Δμs = μs[2:end] .- μs[1:end-1]
ΔkTs = kTs[2:end] .- kTs[1:end-1]
Δμ_ΔkT = Δμs ./ ΔkTs
kTs_centered = (kTs[2:end] .+ kTs[1:end-1]) ./ 2


# Plot the results 
fig = Figure()
xticks = [1, 2, 4, 8]
ax1 = Axis(fig[1,1]; xscale=log10, xticks, xlabel="kT (K)", ylabel="I_Bragg/I_Bragg-max")
ax2 = Axis(fig[1,2]; xscale=log10, xticks, xlabel="kT (K)", ylabel="ΔI_Bragg/ΔkT")
scatter!(ax1, kTs, μs ./ maximum(μs))
scatter!(ax2, kTs_centered, Δμ_ΔkT ./ maximum(μs))


# Inspection of the numerical derivative should show a discontinuity at about 
# T = 3 K. A more careful calculation on a sufficiently large system will show a
# result of T = 3.05 K, which we will use in subsequent scripts.