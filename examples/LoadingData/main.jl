# # Loading Data into Comrade

# The VLBI field does not have a standardized data format, and the EHT uses a
# particular uvfits format similar to the optical interferometry oifits format.
# As a result, we reuse the excellent `eht-imaging` package to load data into `Comrade`.

# Once the data is loaded, we then convert the data into the tabular format `Comrade`
# expects. Note that this may change to a Julia package as the Julia radio
# astronomy group grows.

# To get started, we will load `Comrade` and `Plots` to enable visualizations of the data
import Pkg #hide
__DIR = @__DIR__ #hide
pkg_io = open(joinpath(__DIR, "pkg.log"), "w") #hide
Pkg.activate(__DIR; io=pkg_io) #hide
Pkg.develop(; path=joinpath(__DIR, "..", ".."), io=pkg_io) #hide
Pkg.instantiate(; io=pkg_io) #hide
Pkg.precompile(; io=pkg_io) #hide
close(pkg_io) #hide

ENV["GKSwstype"] = "nul" #hide


using Comrade

using Plots

# We also load Pyehtim since it loads eht-imaging into Julia using PythonCall and exports
# the variable ehtim
using Pyehtim

# To load the data we will use `eht-imaging`. We will use the 2017 public M87 data which can be downloaded from
# [cyverse](https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/EHTC_FirstM87Results_Apr2019)

obseht = ehtim.obsdata.load_uvfits(joinpath(__DIR, "../Data/SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"))
# Now we will average the data over telescope scans. Note that the EHT data has been pre-calibrated so this averaging
# doesn't induce large coherence losses.
obs = Pyehtim.scan_average(obseht)
# !!! warning
#     We use a custom scan-averaging function to ensure that the scan-times are homogenized.
#-
# We can now extract data products that `Comrade` can use
vis    = extract_table(obs, Visibilities()) ## complex visibilites
amp    = extract_table(obs, VisibilityAmplitudes()) ## visibility amplitudes
cphase = extract_table(obs, ClosurePhases(; snrcut=3.0)) ## extract minimal set of closure phases
lcamp  = extract_table(obs, LogClosureAmplitudes(; snrcut=3.0)) ## extract minimal set of log-closure amplitudes

# For polarization we first load the data in the cirular polarization basis
# Additionally, we load the array table at the same time to load the telescope mounts.
obseht = Pyehtim.load_uvfits_and_array(
                joinpath(__DIR, "../Data", "polarized_gaussian_all_corruptions.uvfits"),
                joinpath(__DIR, "..", "Data", "array.txt"),
                polrep="circ"
                        )
obs = Pyehtim.scan_average(obseht)
coh = extract_table(obs, Coherencies())


# !!! warning
#     Always use our `extract_cphase` and `extract_lcamp` functions to find the closures
#     eht-imaging will sometimes incorrectly calculate a non-redundant set of closures.
#-
# We can also recover the array used in the observation using
ac = arrayconfig(vis)
plot(ac) # Plot the baseline coverage

# To plot the data we just call

l = @layout [a b; c d]
pv = plot(vis)
pa = plot(amp)
pcp = plot(cphase)
plc = plot(lcamp)

plot(pv, pa, pcp, plc; layout=l)

# And also the coherency matrices
plot(coh)
