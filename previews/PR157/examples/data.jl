# # Loading Data into Comrade

# The VLBI field does not have a standarized data format, and the EHT uses a
# particular uvfits format that is similar to the optical interferometry oifits format.
# As a result, we reuse the excellent `eht-imaging` package to load data into `Comrade`.

# Once the data is loaded we then convert the data into the tabular format `Comrade`
# expects. Note that in the future this may change to a Julia package as the Julia radio
# astronomy group grows.

# To get started we will load `Comrade` and `Plots` to enable visualizations of the data

using Comrade
using Plots

# Now we load ehtim. This assumes you have a working installation of `eht-imaging`.
# To install eht-imaging see the [ehtim](https://github.com/achael/eht-imaging) github repo.

load_ehtim()

# Now we load the data. We will use the 2017 public M87 data which can be downloaded from
# [cyverse](https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/EHTC_FirstM87Results_Apr2019)

obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "../assets/SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"))
# Add scan and coherently average over them. The eht data has been phase calibrated so that
# this is fine to do.
obs.add_scans()
obs = obs.avg_coherent(0.0, scan_avg=true)

# We can now extract data products that `Comrade` can use
vis = extract_vis(obs) #complex visibilites
amp = extract_amp(obs) # visibility amplitudes
cphase = extract_cphase(obs) # extract minimal set of closure phases
lcamp = extract_lcamp(obs) # extract minimal set of log-closure amplitudes

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
