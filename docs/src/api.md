# Comrade API

## Index

```@index
Pages = ["api.md"]
```

## Model Definitions


### Models

For the description of the model API see [VLBISkyModels](https://ehtjulia.github.io/VLBISkyModels.jl/stable/).




## Data Interface


### Visibility Tables

```@docs
Comrade.ArrayInfo
Comrade.Source
Comrade.stationindex
Comrade.VisTable
Comrade.Correlation
Comrade.Correlations
Comrade.measurement(::Comrade.VisTable)
Comrade.weight
Comrade.noise(::Comrade.VisTable)
Comrade.domain(::Comrade.VisTable)
Comrade.info(::Comrade.VisTable)
Comrade.stationgeometry(::Comrade.VisTable)
Comrade.baselines
Comrade.stations(::Comrade.VisTable)
Comrade.times
Comrade.frequencies
Comrade.beamsize(::Comrade.VisTable)
Comrade.amplitudes
Comrade.flag
Comrade.drop
Comrade.add_fractional_noise
Comrade.load_array_txt
Comrade.reset_mounts
Comrade.datarows
Comrade.ecef_to_geodetic
Comrade.elevation_parallactic
```

### Averaging

```@docs
Comrade.average
Comrade.ByScan
Comrade.ChannelBins
Comrade.ByChannel
Comrade.averaging_axis
```

### Closures

```@docs
Comrade.ClosureTable
Comrade.closure_phases(::Comrade.VisTable)
Comrade.logclosure_amplitudes
Comrade.noise(::Comrade.ClosureTable)
Comrade.ClosurePhases
Comrade.LogClosureAmplitudes
```



## VLBI Modeling

### Sky Models

```@docs
Comrade.@sky
Comrade.AbstractSkyModel
Comrade.SkyModel
Comrade.FixedSkyModel
Comrade.MultiSkyModel
Comrade.idealmaps
Comrade.skymodel(::Comrade.AbstractVLBIPosterior, ::Any)
```

### Instrument Models

```@docs
Comrade.@instrument
Comrade.@jones
Comrade.CalTable
Comrade.caltable(::Comrade.StationArray)
Comrade.stationparams
Comrade.stations(::Comrade.CalTable)
Comrade.IIDStationPrior
Comrade.GaussMarkovStationPrior
Comrade.AbstractChainParameterization
Comrade.Centered
Comrade.NonCentered
Comrade.AngleEmbedded
Comrade.AbstractGaussMarkovProcess
Comrade.OrnsteinUhlenbeck
Comrade.BrownianMotion
Comrade.WrappedBrownian
Comrade.WrappedOrnsteinUhlenbeck
Comrade.VonMisesProcess
Comrade.stationary_moments
Comrade.transition_moments
Comrade.isstationary
Comrade.AbstractInitialPrior
Comrade.StationaryInit
Comrade.FixedInit
Comrade.GaussianInit
Comrade.UniformInit
Comrade.initial_moments
Comrade.marginal_moments
Comrade.hyperprior
Comrade.materialize
Comrade.ArrayPrior
Comrade.Segmentation
Comrade.IntegSeg
Comrade.ScanSeg
Comrade.TrackSeg
Comrade.EpochSeg
Comrade.timestamps
Comrade.SingleReference
Comrade.MultiReference
Comrade.SensitiveReference
Comrade.StationPhase
Comrade.ReportGauge
Comrade.PinGauge
Comrade.gauge_pins
Comrade.default_gauge_preference
Comrade.AbstractJonesMatrix
Comrade.Jones
Comrade.gains
Comrade.leakage
Comrade.rotation
Comrade.JonesProduct
Comrade.PerStation
Comrade.AbstractInstrumentModel
Comrade.IdealInstrumentModel
Comrade.InstrumentModel
Comrade.StationArray
Comrade.StationLookup
Comrade.forward_jones
Comrade.selfcal_phase
```


### Posterior Constructions

```@docs
Comrade.AbstractVLBIPosterior
Comrade.logprior
Comrade.loglikelihood
Comrade.dataproducts
Comrade.skymodel
Comrade.instrumentmodel(::Comrade.AbstractVLBIPosterior)
Comrade.instrumentmodel(::Comrade.AbstractVLBIPosterior, ::Any)
Comrade.prior_sample
Comrade.likelihood
Comrade.VLBIPosterior
Comrade.simulate_observation
Comrade.residuals
Comrade.TransformedVLBIPosterior
Comrade.ImgNormalData
Comrade.datadomain
Comrade.VisDomain
Comrade.ImageDomain
Comrade.makelikelihood
Comrade.transform(::Comrade.TransformedVLBIPosterior, ::Any)
Comrade.inverse(::Comrade.TransformedVLBIPosterior, ::Any)
Comrade.ascube(::Comrade.VLBIPosterior)
Comrade.asflat(::Comrade.VLBIPosterior)
Comrade.transport_to(::Comrade.VLBIPosterior, ::Comrade.PT.AbstractStdDist)
```

### Inference
```@docs
Comrade.estimate_mode
Comrade.MemoryStore
Comrade.DiskStore
Comrade.load_samples
Comrade.PosteriorSamples
Comrade.postsamples
Comrade.samplerstats
Comrade.samplerinfo
Comrade.resample_equal
Comrade.residualplot
Comrade.chi2
Comrade.ndata
Comrade.ReactantNUTS
Comrade.default_disk_callback
```

### Moves

```@docs
Comrade.AbstractMove
Comrade.MoveSet
Comrade.NormalStep
Comrade.ComponentwiseNormalStep
Comrade.SymmetricDiscreteStep
Comrade.CompensatedMove
Comrade.ChainHyperMove
Comrade.PhaseSheetMove
Comrade.chain_hyper_moves
Comrade.flux_gain_move
Comrade.CoordinateView
Comrade.check_move
Comrade.check_moveset
Comrade.move_summary
Comrade.move_seconds
```

The protocol a move implements, and the helpers for writing and checking one:

```@docs
Comrade.propose
Comrade.step_distribution
Comrade.move_name
Comrade.likelihood_invariant
Comrade.context_data
Comrade.traceable
Comrade.move_context
Comrade.update_context!
Comrade.step_label
Comrade.component_coords
Comrade.component_label
Comrade.node
Comrade.coords
Comrade.value
Comrade.latent
```

## Misc

```@docs
Comrade.station_tuple
Comrade.dirty_image
Comrade.dirty_beam
Comrade.beamsize
Comrade.apply_fluctuations
Comrade.apply_fluctuations!
Comrade.UnitFluxMap
Comrade.corr_image_prior
Comrade.rmap
```

## Internal (Not Public API)

```@docs
Comrade.ObservedSkyModel
Comrade.forward_model_map
Comrade.fixed_indices

```


## Plotting


**Warning**
A user must first load `Makie` or a `Makie` backend, e.g., `CairoMakie` to use this functionality

```@docs
Comrade.fieldplot
Comrade.datafield
Comrade.DataField
Comrade.caltableplot
```