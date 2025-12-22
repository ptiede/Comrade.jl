# Comrade API

## Index

```@index
Pages = ["api.md"]
```

```@docs
Comrade.Comrade
```

## Model Definitions


### Models

For the description of the model API see [VLBISkyModels](https://ehtjulia.github.io/VLBISkyModels.jl/stable/).




## Data Interface


### Data Tables

```@docs
Comrade.AbstractVLBITable
Comrade.datatable
Comrade.AbstractArrayConfiguration
Comrade.EHTArrayBaselineDatum
Comrade.EHTArrayConfiguration
Comrade.ClosureConfig
Comrade.sites(::Comrade.AbstractArrayConfiguration)
Comrade.domain(::Comrade.AbstractArrayConfiguration)
Comrade.beamsize(::Comrade.AbstractArrayConfiguration)
Comrade.logclosure_amplitudes
Comrade.closure_phases
Comrade.AbstractObservationTable
Comrade.EHTObservationTable
Comrade.domain(::Comrade.AbstractObservationTable)
Comrade.arrayconfig(::Comrade.AbstractObservationTable)
Comrade.beamsize(::Comrade.AbstractObservationTable)
Comrade.sites(::Comrade.AbstractObservationTable)
Comrade.TimeTable
Comrade.Scan
Comrade.timetable
Comrade.flag
Base.filter(::Any, ::Comrade.EHTObservationTable)
Comrade.select_baseline
Comrade.add_fractional_noise
```

### Datums

```@docs
Comrade.AbstractVisibilityDatum
Comrade.EHTCoherencyDatum
Comrade.EHTVisibilityDatum
Comrade.EHTVisibilityAmplitudeDatum
Comrade.EHTLogClosureAmplitudeDatum
Comrade.EHTClosurePhaseDatum
Comrade.triangle
Comrade.baseline
Comrade.quadrangle
```

### Data Products

```@docs
Comrade.extract_table
Comrade.Visibilities
Comrade.VisibilityAmplitudes
Comrade.ClosurePhases
Comrade.LogClosureAmplitudes
Comrade.Coherencies
```



## VLBI Modeling

### Sky Models

```@docs
Comrade.AbstractSkyModel
Comrade.SkyModel
Comrade.FixedSkyModel
Comrade.MultiSkyModel
Comrade.idealmaps
Comrade.skymodel(::Comrade.AbstractVLBIPosterior, ::Any)
```

### Instrument Models

```@docs
Comrade.CalTable
Comrade.caltable(::Comrade.SiteArray)
Comrade.sites(::Comrade.CalTable)
Comrade.IIDSitePrior
Comrade.ArrayPrior
Comrade.Segmentation
Comrade.IntegSeg
Comrade.ScanSeg
Comrade.TrackSeg
Comrade.timestamps
Comrade.SingleReference
Comrade.SEFDReference
Comrade.SingleStokesGain
Comrade.JonesG
Comrade.JonesD
Comrade.JonesR
Comrade.JonesF
Comrade.GenericJones
Comrade.JonesSandwich
Comrade.AbstractInstrumentModel
Comrade.IdealInstrumentModel
Comrade.InstrumentModel
Comrade.SiteArray
Comrade.SiteLookup
Comrade.forward_jones
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
Comrade.forward_model
Comrade.prior_sample
Comrade.likelihood
Comrade.VLBIPosterior
Comrade.simulate_observation
Comrade.residuals
Comrade.TransformedVLBIPosterior
HypercubeTransform.transform(::Comrade.TransformedVLBIPosterior, ::Any)
HypercubeTransform.inverse(::Comrade.TransformedVLBIPosterior, ::Any)
HypercubeTransform.ascube(::Comrade.VLBIPosterior)
HypercubeTransform.asflat(::Comrade.VLBIPosterior)
```

### Inference
```@docs
Comrade.comrade_opt
Comrade.MemoryStore
Comrade.DiskStore
Comrade.load_samples
Comrade.PosteriorSamples
Comrade.postsamples
Comrade.samplerstats
Comrade.samplerinfo
Comrade.resample_equal
Comrade.residual
Comrade.residual_data
Comrade.chi2
```

## Misc

```@docs
Comrade.site_tuple
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
Comrade.build_datum
Comrade.ObservedSkyModel
```

### eht-imaging interface (Internal)

```@docs
Comrade.extract_amp
Comrade.extract_cphase
Comrade.extract_lcamp
Comrade.extract_vis
Comrade.extract_coherency
```

## Plotting


**Warning**
A user must first load `Makie` or a `Makie` backend, e.g., `CairoMakie` to use this functionality

```@docs
Comrade.plotfields
Comrade.plotfields!
Comrade.axisfields
Comrade.plotcaltable
Comrade.baselineplot
Comrade.baselineplot!
ComradeMakieExt.BaselinePlot
```