using LinearAlgebra

# Physical basis for SU(3). Needed for calculating the generalized classical sum rule.
function observable_matrices(; N=3)
    Sx, Sy, Sz = Sunny.spin_matrices(; N)
    Os = [
        Sx,
        Sy,
        Sz,
        -(Sx*Sz + Sz*Sx),
        -(Sy*Sz + Sz*Sy),
        Sx^2 - Sy^2,
        Sx*Sy + Sy*Sx,
        √3 * Sz^2 - I*2/√3,
    ] 
    return Os
end


# Calculate the the integral of S(𝐪, ω) over both 𝐪 and ω.
function total_spectral_weight(sc::SampledCorrelations; kT = Inf)
    qs = available_wave_vectors(sc)
    formula = intensity_formula(sc, :trace; kT)
    is = intensities_interpolated(sc, qs, formula; interpolation=:round, negative_energies=true)
    return sum(is)
end


# Renormalize kets and adjust κs in system so renormalization is applied
# when running the dynamics.
function renormalize_system!(sys, coherents, κ)
    sys.κs .= κ
    for site in Sunny.eachsite(sys)
        set_coherent!(sys, coherents[site], site)
    end
    return nothing
end
renormalize_system!(sys, κ) = renormalize_system!(sys, sys.coherents, κ)


# Estimate S(𝐪,ω) at temperature kT and evaluate the sum rule using both the
# classical-to-quantum correspondence factor and moment renormalization. For the
# publication results, the intensities at the ordering wave vector were
# integrated out prior to applying the classical-to-quantum correspondence
# factor whenever T < T_N. In other words, the classical-to-quantum rescaling
# was only applied to the inelastic response.
function estimate_sum(sys::System{N}, κ, kT, sim_params; observables=nothing) where N
    (; nω, ωmax, Δt, Δt_therm, dur_therm, dur_decorr, nsamples, λ) = sim_params

    sc = dynamical_correlations(sys; Δt, nω, ωmax, observables)
    saved_coherents = copy(sys.coherents)
    ndecorr = round(Int, dur_decorr/Δt_therm)
    ntherm = round(Int, dur_therm/Δt_therm)

    # Thermalize the system
    langevin = Langevin(Δt_therm; λ, kT = kT*Sunny.meV_per_K)
    for _ in 1:ntherm
        step!(sys, langevin)
    end

    # Collect samples
    for _ in 1:nsamples

        # Get a decorrelated equilibrium sample 
        for _ in 1:ndecorr
            step!(sys, langevin)
        end

        # Save the sampled state
        saved_coherents .= sys.coherents

        # Renormalize kets before running dynamics to collect a trajectory and associated correlations.
        renormalize_system!(sys, sys.coherents, κ)
        add_sample!(sc, sys)

        # Un-renormalize kets and return to previous state before continuing
        # Langevin sampling process in next iteration.
        renormalize_system!(sys, saved_coherents, 1.0)
    end

    # Return spectral weight per site.
    return total_spectral_weight(sc; kT = kT*Sunny.meV_per_K) / length(Sunny.eachsite(sys))
end