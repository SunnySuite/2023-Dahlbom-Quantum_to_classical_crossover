using DrWatson
@quickactivate "FiniteTemperatureFeI2"

using Sunny, GLMakie
include(srcdir("instrument_corrections.jl"))

# # Examining intensity data
#
# Here intenity information is extracted from the structure factors calculated
# in the previous script.
data_renormalized = load(datadir("dssf", "renormalized.jld2"))
data_unrenormalized = load(datadir("dssf", "unrenormalized.jld2"))

# ## "Spaghetti" plots
# Define a path through reciprocal space. 
cryst = data_renormalized["cryst"]
points = [[0,   0, 0],  # List of wave vectors that define a path
          [1,   0, 0],
          [0,   1, 0],
          [1/2, 0, 0],
          [0,   1, 0],
          [0,   0, 0]] 
density = 60
path, xticks = reciprocal_space_path(cryst, points, density)

# Set up figure for spaghetti plots,
fig = Figure()
ax1 = Axis(fig[1,1]; xticks, xlabel = "[H,K,L]", ylabel="Energy (meV)", title="Unrenormalized")
ax2 = Axis(fig[1,2]; xticks, xlabel = "[H,K,L]", ylabel="Energy (meV)", title="Renormalized")

# and plot both the unrenormalized results,
sc = data_unrenormalized["sc"]
energies = available_energies(sc)
formula = intensity_formula(sc, :perp; kT=kT_meV, formfactors = [FormFactor("Fe2"; g_lande=3/2)])

is = intensities_interpolated(sc, path, formula)
is = convolve_sequoia(is, energies)  # Add instrumental broadening for Sequoia instrument at SNS
heatmap!(ax1, 1:size(is, 1), energies, is; colorrange=(0.0, 0.5))

# and the renormalized results.
sc_rn = data_renormalized["sc"]

is = intensities_interpolated(sc_rn, path, formula)
is = convolve_sequoia(is, energies)  # Add instrumental broadening
heatmap!(ax2, 1:size(is, 1), energies, is; colorrange=(0.0, 0.5))
fig


# ## Single-𝐪 cuts
# 
# The approach is similar to above, except we will look at one wave vector at a
# time.
qs = [[1/2, 0, 0], [3/4, 0, 0], [1, 0, 0]]        # Wave vectors to examine 
σ = σ_from_FWHM(0.47)                             # Broadening width for HB3 triple-axis spectrometer at HIFER
gauss(x, x₀) = 1/(σ*√2π)*exp(-(x-x₀)^2 / (2σ^2))  # Broadening function

cuts = intensities_interpolated(sc, qs, formula)
cuts = broaden_energy(sc, cuts, gauss)

cuts_rn = intensities_interpolated(sc_rn, qs, formula)
cuts_rn = broaden_energy(sc, cuts_rn, gauss)

# Plot comparisons at individual wave vectors.
fig = Figure()

ax1 = Axis(fig[1,1]; ylabel = "Intensity (a.u.)", xlabel = "Energy (meV)", title="Q=(1/2, 0, 0)")
ax2 = Axis(fig[1,2]; xlabel = "Energy (meV)", title="Q=(3/4, 0, 0)")
ax3 = Axis(fig[1,3]; xlabel = "Energy (meV)", title="Q=(1, 0, 0)")
color = :blue

lines!(ax1, energies, cuts[1,:];    linestyle=:dash, color)
lines!(ax1, energies, cuts_rn[1,:]; color)
lines!(ax2, energies, cuts[2,:];    linestyle=:dash, color)
lines!(ax2, energies, cuts_rn[2,:]; color)
ln = lines!(ax3, energies, cuts[3,:];    linestyle=:dash, color)
ln_rn = lines!(ax3, energies, cuts_rn[3,:]; color)

Legend(fig[1,4], [ln, ln_rn], ["Unrenromalized", "Renormalized"])
fig