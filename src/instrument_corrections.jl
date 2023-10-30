σ_from_FWHM(FWHM) = FWHM / (2*√2log(2))

function convolve_sequoia(slice, ωs)
    @assert length(ωs) == size(slice)[end] "Number of enegy values does not match slice dimension"

    # Load sequoia resolution data
    data = CSV.read(datadir("instrument_data", "sequoia_FWHM_per_energy.csv"), DataFrame)
    Es, FWHMs = data.energy, data.FWHM
    σs = σ_from_FWHM.(FWHMs)
    println(σs)
    σ_seq = LinearInterpolation(Es, σs)

    slice′ = zero(slice)
    for (i, ω) ∈ enumerate(ωs)
        if ω > 11.3
            error("No resolution data available above 11.3 meV")
        end
        σ = σ_seq[ω]
        gauss = @. 1/(σ*√(2π))*exp(-(ωs - ω)^2 / (2σ^2))
        for j ∈ 1:size(slice)[1]
            slice′[j,i] = slice[j,:] ⋅ gauss
        end
    end

    return slice′
end