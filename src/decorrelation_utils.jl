# Run system at fixed temperature for specified duration.
function thermalize!(sys, kT, Δt, dur; λ = 0.1)
    nsteps = round(Int, dur/Δt)
    integrator = Langevin(Δt; λ, kT)
    for _ in 1:nsteps
        step!(sys, integrator)
    end
end

# Lookup decorrelation time and run dynamics for this duration.
function decorrelate!(sys, kT, Δt; λ=0.1)
    integrator = Langevin(Δt; kT, λ)
    nsteps = round(Int, decorrelation_time(kT)/Δt)
    for _ in 1:nsteps
        step!(sys, integrator)
    end
end

# Estimated by examining the decorrelation time of the energy
# time series generated at different temperatures.
function decorrelation_time(kT)
    Δt = 0.004
    @assert kT > 0
    if 0 < kT < 1.13
        3000 * Δt
    elseif kT < 1.83
        1000.0 * Δt
    elseif kT < 3.79
        600.0 * Δt
    elseif kT < 4.83
        400.0 * Δt
    elseif kT < 6.16 
        350.0 * Δt
    elseif kT < 7.85
        300.0 * Δt
    else
        200.0 * Δt
    end
end