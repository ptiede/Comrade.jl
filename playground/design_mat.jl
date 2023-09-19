using Comrade
using Plots
using Pyehtim

# include("/Users/urirolls/Comrade.jl/src/calibration/jones.jl")

obseht = ehtim.obsdata.load_uvfits(joinpath(dirname(pathof(Comrade)), "..", "examples", "SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"))
obs = Pyehtim.scan_average(obseht)
vis = extract_table(obs, ComplexVisibilities()) #complex visibilites

mat = jonescache(vis, ScanSeg())


# function timetable(obs::EHTObservation)
#     obs = vis
#     times = obs[:T]
#     sites = stations(obs) 
#     bs = obs[:baseline]
#     scantimes = Vector{Vector{Float64}}()
#     for site in sites
#         mask = [site in n_bs for n_bs in bs]
#         filtered_values = [v for (v, m) in zip(times, mask) if m]
#         site_times = unique(filtered_values)
#         push!(scantimes, site_times)
#         TimeTable()
#     end

# struct TimeTable{O, S,T}
#     obs:: O
#     sites:: S
#     times:: T
# end