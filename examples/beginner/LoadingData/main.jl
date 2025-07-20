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
# As a result, we reuse the excellent `eht-imaging` package to load data into `Comrade`.

# Once the data is loaded, we then convert the data into the tabular format `Comrade`
# expects. Note that this may change to a Julia package as the Julia radio
# astronomy group grows.

# To get started, we will load `Comrade` and `Plots` to enable visualizations of the data
using Comrade
using CairoMakie

# We also load Pyehtim since it loads eht-imaging into Julia using PythonCall and exports
# the variable ehtim
using Pyehtim

# To load the data we will use `eht-imaging`. We will use the 2017 public M87 data which can be downloaded from
# [cyverse](https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/EHTC_FirstM87Results_Apr2019)

obseht = ehtim.obsdata.load_uvfits(joinpath(__DIR, "..", "..", "Data", "SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"))
# Now we will average the data over telescope scans. Note that the EHT data has been pre-calibrated so this averaging
# doesn't induce large coherence losses.
obs = Pyehtim.scan_average(obseht)
# !!! warning
#     We use a custom scan-averaging function to ensure that the scan-times are homogenized.
#-
# We can now extract data products that `Comrade` can use
vis = extract_table(obs, Visibilities()) ## complex visibilites
amp = extract_table(obs, VisibilityAmplitudes()) ## visibility amplitudes
cphase = extract_table(obs, ClosurePhases(; snrcut = 3.0)) ## extract minimal set of closure phases
lcamp = extract_table(obs, LogClosureAmplitudes(; snrcut = 3.0)) ## extract minimal set of log-closure amplitudes

# For polarization we first load the data in the cirular polarization basis
# Additionally, we load the array table at the same time to load the telescope mounts.
obseht = Pyehtim.load_uvfits_and_array(
    joinpath(__DIR, "..", "..", "Data", "polarized_gaussian_all_corruptions.uvfits"),
    joinpath(__DIR, "..", "..", "Data", "array.txt"),
    polrep = "circ"
)
obs = Pyehtim.scan_average(obseht)
coh = extract_table(obs, Coherencies())


# !!! warning
#     Always use our `extract_cphase` and `extract_lcamp` functions to find the closures
#     eht-imaging will sometimes incorrectly calculate a non-redundant set of closures.
#-
# We can also recover the array used in the observation using
using DisplayAs
plotfields(coh, U, V, axis_kwargs = (xreversed = true,)) |> DisplayAs.PNG |> DisplayAs.Text # Plot the baseline coverage

# As of Comrade 0.11.7 Makie is the preferred plotting tool. For plotting data there are two
# classes of functions:
#  - `baselineplot` which gives complete control of plotting
#  - `plotfields, axisfields` which are more automated and limited but will automatically add
#     labels, legends, titles etc.
fig = Figure(; size = (800, 600))
plotfields!(fig[1, 1], vis, uvdist, measurement)
plotfields!(fig[1, 2], amp, uvdist, measurement)
plotfields!(fig[2, 1], cphase, uvdist, measurement)
plotfields!(fig[2, 2], lcamp, uvdist, measurement)
fig |> DisplayAs.PNG |> DisplayAs.Text

# And also the coherency matrices. Since the data products are a matrix we need to plot each one separately.
fig = Figure(; size = (800, 600))
plotfields!(fig[1, 1], coh, uvdist, x -> measwnoise(x)[1, 1], axis_kwargs = (ylabel = "RR", xlabel = "uv distance (Gλ)"))
plotfields!(fig[2, 1], coh, uvdist, x -> measwnoise(x)[2, 1], axis_kwargs = (ylabel = "LR", xlabel = "uv distance (Gλ)"))
plotfields!(fig[1, 2], coh, uvdist, x -> measwnoise(x)[1, 2], axis_kwargs = (ylabel = "RL", xlabel = "uv distance (Gλ)"))
plotfields!(fig[2, 2], coh, uvdist, x -> measwnoise(x)[2, 2], axis_kwargs = (ylabel = "LL", xlabel = "uv distance (Gλ)"))
fig

# You can also plot a single baseline
fig, ax = plotfields(coh, (:AA, :LM), Ti, x -> abs(measwnoise(x)[1, 1]), axis_kwargs = (; ylabel = "|RR|"))
ax2 = plotfields!(fig[1, 2], coh, (:LM, :AZ), Ti, x -> abs(measwnoise(x)[1, 1]), axis_kwargs = (; ylabel = "|RR|"))
fig

# Finally, we provide a more low-level plotting function `baselineplot` which allows you to plot
# any field against any other field. This is what `plotfields` calls under the hood. However, it
# does not automatically add labels, legends, titles etc, but can add multiple baselines to the same plot.
fig, ax = baselineplot(coh, (:AA, :LM), Ti, x -> abs(measwnoise(x)[1, 1]), label = "AA-LM")
baselineplot!(ax, coh, (:LM, :AZ), Ti, x -> abs(measwnoise(x)[1, 1]), label = "LM-AZ")
ax.ylabel = "|RR|"
ax.xlabel = "Time (UTC)"
axislegend(ax)
fig
