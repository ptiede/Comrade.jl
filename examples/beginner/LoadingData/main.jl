import Pkg #hide
__DIR = @__DIR__ #hide
pkg_io = open(joinpath(__DIR, "pkg.log"), "w") #hide
Pkg.activate(__DIR; io = pkg_io) #hide
Pkg.develop(; path = joinpath(__DIR, "..", "..", ".."), io = pkg_io) #hide
Pkg.instantiate(; io = pkg_io) #hide
Pkg.precompile(; io = pkg_io) #hide
close(pkg_io) #hide

# # Loading Data into Comrade

# The VLBI field does not have a standardized data format, and the EHT uses a
# particular uvfits format similar to the optical interferometry oifits format.
# In Comrade we read uvfits with the pure-Julia [`VLBIFiles.jl`](https://github.com/JuliaAstro/VLBIFiles.jl)
# package, which avoids any Python dependency.

# Once the data is loaded, we then convert the data into the tabular format `Comrade`
# expects.

# To get started, we will load `Comrade` and `CairoMakie` to enable visualizations of the data.
using Comrade
using CairoMakie

# We also load `VLBIFiles` to read the uvfits file.
using VLBIFiles

# We will use the 2017 public M87 data which can be downloaded from
# [cyverse](https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/EHTC_FirstM87Results_Apr2019)
uvd = VLBIFiles.load(
    VLBIFiles.UVData,
    joinpath(__DIR, "..", "..", "Data", "SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits")
)

# `VisTable(uvd, IPol())` is the table of Stokes I visibilities, and `average(t, ByScan())`
# averages it over telescope scans. Note that the EHT data has been pre-calibrated so this
# averaging doesn't induce large coherence losses.
vis = average(VisTable(uvd, IPol()), ByScan())   # complex visibilities

# The other data products are formed from the averaged visibilities.
amp = amplitudes(vis)                    # visibility amplitudes
cphase = closure_phases(vis)             # minimal closure phases
lcamp = log_closure_amplitudes(vis)      # minimal log-closure amplitudes

# For polarization we load the file the same way; circular polarization is detected from the
# antenna feed types in the FITS table. `VisTable(uvdp, Correlations())` holds all four feed
# products of each baseline as a 2×2 coherency matrix, and `reset_mounts!` sets the telescope
# mount information from the antenna text file.
uvdp = VLBIFiles.load(
    VLBIFiles.UVData,
    joinpath(__DIR, "..", "..", "Data", "polarized_gaussian_all_corruptions.uvfits")
)
coh = average(VisTable(uvdp, Correlations()), ByScan())
reset_mounts!(coh, joinpath(__DIR, "..", "..", "Data", "array.txt"))

#-
# We can also recover the array used in the observation using
using DisplayAs
plotfields(coh, U, V, axis_kwargs = (xreversed = true,)) |> DisplayAs.PNG |> DisplayAs.Text # Plot the baseline coverage

# As of Comrade 0.11.7 Makie is the preferred plotting tool. For plotting data there are two
# classes of functions:
#  - `baselineplot` which gives complete control of plotting
#  - `plotfields, axisfields` which are more automated and limited but will automatically add
#     labels, legends, titles etc.
#
# For the closures, `:uvdist` is the longest baseline in each closure.
fig = Figure(; size = (800, 600))
plotfields!(fig[1, 1], vis, :uvdist, :measurement)
plotfields!(fig[1, 2], amp, :uvdist, :measurement)
plotfields!(fig[2, 1], cphase, :uvdist, :measurement)
plotfields!(fig[2, 2], lcamp, :uvdist, :measurement)
fig |> DisplayAs.PNG |> DisplayAs.Text

# And also the coherency matrices. Since the data products are a matrix we need to plot each one
# separately. Unmeasured feed products are `NaN`.
fig = Figure(; size = (800, 600))
plotfields!(fig[1, 1], coh, :uvdist, x -> x.measurement[1, 1], axis_kwargs = (ylabel = "RR", xlabel = "uv distance (λ)"))
plotfields!(fig[2, 1], coh, :uvdist, x -> x.measurement[2, 1], axis_kwargs = (ylabel = "LR", xlabel = "uv distance (λ)"))
plotfields!(fig[1, 2], coh, :uvdist, x -> x.measurement[1, 2], axis_kwargs = (ylabel = "RL", xlabel = "uv distance (λ)"))
plotfields!(fig[2, 2], coh, :uvdist, x -> x.measurement[2, 2], axis_kwargs = (ylabel = "LL", xlabel = "uv distance (λ)"))
fig

# You can also plot a single baseline
fig, ax = plotfields(coh, (:AA, :LM), Ti, x -> abs(x.measurement[1, 1]), axis_kwargs = (; ylabel = "|RR|"))
ax2 = plotfields!(fig[1, 2], coh, (:LM, :AZ), Ti, x -> abs(x.measurement[1, 1]), axis_kwargs = (; ylabel = "|RR|"))
fig

# Finally, we provide a more low-level plotting function `baselineplot` which allows you to plot
# any field against any other field. This is what `plotfields` calls under the hood. However, it
# does not automatically add labels, legends, titles etc, but can add multiple baselines to the same plot.
fig, ax = baselineplot(coh, (:AA, :LM), Ti, x -> abs(x.measurement[1, 1]), label = "AA-LM")
baselineplot!(ax, coh, (:LM, :AZ), Ti, x -> abs(x.measurement[1, 1]), label = "LM-AZ")
ax.ylabel = "|RR|"
ax.xlabel = "Time (UTC)"
axislegend(ax)
fig
