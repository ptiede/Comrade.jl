# Stoked API

## Index

```@index
Pages = ["api.md"]
```

```@docs
Stoked.Stoked
```

## Model Definitions


### Models

For the description of the model API see [VLBISkyModels](https://ehtjulia.github.io/VLBISkyModels.jl/stable/).




## Data Interface


### Data Tables

```@docs
Stoked.AbstractVLBITable
Stoked.datatable
Stoked.AbstractArrayConfiguration
Stoked.EHTArrayBaselineDatum
Stoked.EHTArrayConfiguration
Stoked.ClosureConfig
Stoked.sites(::Stoked.AbstractArrayConfiguration)
Stoked.domain(::Stoked.AbstractArrayConfiguration)
Stoked.beamsize(::Stoked.AbstractArrayConfiguration)
Stoked.logclosure_amplitudes
Stoked.closure_phases
Stoked.AbstractObservationTable
Stoked.EHTObservationTable
Stoked.domain(::Stoked.AbstractObservationTable)
Stoked.arrayconfig(::Stoked.AbstractObservationTable)
Stoked.beamsize(::Stoked.AbstractObservationTable)
Stoked.sites(::Stoked.AbstractObservationTable)
Stoked.TimeTable
Stoked.Scan
Stoked.timetable
Stoked.flag
Base.filter(::Any, ::Stoked.EHTObservationTable)
Stoked.select_baseline
Stoked.add_fractional_noise
```

### Datums

```@docs
Stoked.AbstractVisibilityDatum
Stoked.EHTCoherencyDatum
Stoked.EHTVisibilityDatum
Stoked.EHTVisibilityAmplitudeDatum
Stoked.EHTLogClosureAmplitudeDatum
Stoked.EHTClosurePhaseDatum
Stoked.triangle
Stoked.baseline
Stoked.quadrangle
```

### Data Products

```@docs
Stoked.extract_table
Stoked.Visibilities
Stoked.VisibilityAmplitudes
Stoked.ClosurePhases
Stoked.LogClosureAmplitudes
Stoked.Coherencies
```



## VLBI Modeling

### Sky Models

```@docs
Stoked.AbstractSkyModel
Stoked.SkyModel
Stoked.FixedSkyModel
Stoked.idealvisibilities
Stoked.skymodel(::Stoked.AbstractVLBIPosterior, ::Any)
```

### Instrument Models

```@docs
Stoked.CalTable
Stoked.caltable(::Stoked.SiteArray)
Stoked.sites(::Stoked.CalTable)
Stoked.IIDSitePrior
Stoked.ArrayPrior
Stoked.Segmentation
Stoked.IntegSeg
Stoked.ScanSeg
Stoked.TrackSeg
Stoked.timestamps
Stoked.SingleReference
Stoked.SEFDReference
Stoked.SingleStokesGain
Stoked.JonesG
Stoked.JonesD
Stoked.JonesR
Stoked.JonesF
Stoked.GenericJones
Stoked.JonesSandwich
Stoked.AbstractInstrumentModel
Stoked.IdealInstrumentModel
Stoked.InstrumentModel
Stoked.SiteArray
Stoked.SiteLookup
Stoked.forward_jones
```


### Posterior Constructions

```@docs
Stoked.AbstractVLBIPosterior
Stoked.logprior
Stoked.loglikelihood
Stoked.dataproducts
Stoked.skymodel
Stoked.instrumentmodel(::Stoked.AbstractVLBIPosterior)
Stoked.instrumentmodel(::Stoked.AbstractVLBIPosterior, ::Any)
Stoked.forward_model
Stoked.prior_sample
Stoked.likelihood
Stoked.VLBIPosterior
Stoked.simulate_observation
Stoked.residuals
Stoked.TransformedVLBIPosterior
HypercubeTransform.transform(::Stoked.TransformedVLBIPosterior, ::Any)
HypercubeTransform.inverse(::Stoked.TransformedVLBIPosterior, ::Any)
HypercubeTransform.ascube(::Stoked.VLBIPosterior)
HypercubeTransform.asflat(::Stoked.VLBIPosterior)
```

### Inference
```@docs
Stoked.comrade_opt
Stoked.MemoryStore
Stoked.DiskStore
Stoked.load_samples
Stoked.PosteriorSamples
Stoked.postsamples
Stoked.samplerstats
Stoked.samplerinfo
Stoked.resample_equal
Stoked.residual
Stoked.residual_data
Stoked.chi2
```

## Misc

```@docs
Stoked.site_tuple
Stoked.dirty_image
Stoked.dirty_beam
Stoked.beamsize
Stoked.apply_fluctuations
Stoked.corr_image_prior
Stoked.rmap
```

## Internal (Not Public API)

```@docs
Stoked.build_datum
Stoked.ObservedSkyModel
```

### eht-imaging interface (Internal)

```@docs
Stoked.extract_amp
Stoked.extract_cphase
Stoked.extract_lcamp
Stoked.extract_vis
Stoked.extract_coherency
```

## Plotting


**Warning**
A user must first load `Makie` or a `Makie` backend, e.g., `CairoMakie` to use this functionality

```@docs
Stoked.plotfields
Stoked.axisfields
Stoked.plotcaltable
Stoked.baselineplot
```