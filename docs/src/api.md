# Comrade API



## Contents

```@contents
Pages = ["api.md"]
```

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
Comrade.datatable(::Comrade.AbstractVLBITable)
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
Comrade.datatable(::Comrade.AbstractObservationTable)
Comrade.domain(::Comrade.AbstractObservationTable)
Comrade.arrayconfig(::Comrade.AbstractObservationTable)
Comrade.beamsize(::Comrade.AbstractObservationTable)
Comrade.sites(::Comrade.AbstractObservationTable)
Comrade.TimeTable
Comrade.Scan
Comrade.timetable
```

### Datums

```@docs
Comrade.AbstractVisibilityDatum
Comrade.EHTCoherencyDatum
Comrade.EHTVisibilityDatum
Comrade.EHTVisibilityAmplitudeDatum
Comrade.EHTLogClosureAmplitudeDatum
Comrade.EHTClosurePhaseDatum
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
```

### Instrument Models

```@docs
Comrade.CalTable
Comrade.caltable(::Comrade.SiteArray)
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
Comrade.JonesSandwich
Comrade.AbstractInstrumentModel
Comrade.IdealInstrumentModel
Comrade.InstrumentModel
Comrade.SiteArray
Comrade.SiteLookup
```


### Posterior Constructions

```@docs
Comrade.AbstractVLBIPosterior
Comrade.logprior
Comrade.loglikelihood
Comrade.dataproducts
Comrade.skymodel(::Comrade.AbstractVLBIPosterior)
Comrade.instrumentmodel(::Comrade.AbstractVLBIPosterior)
Comrade.forward_model
Comrade.prior_sample
Comrade.likelihood
Comrade.VLBIPosterior
Comrade.simulate_observation
Comrade.residuals
Comrade.TransformedVLBIPosterior
HypercubeTransform.transform(::Comrade.AbstractVLBIPosterior, ::Any)
HypercubeTransform.inverse(::Comrade.AbstractVLBIPosterior, ::Any)
HypercubeTransform.ascube(::Comrade.AbstractVLBIPosterior)
HypercubeTransform.asflat(::Comrade.AbstractVLBIPosterior)
```

### Inference
```@docs
Comrade.comrade_opt
Dynesty.dysample(::Comrade.VLBIPosterior)
AbstractMCMC.sample(rng::AbstractRNG, ::Comrade.VLBIPosterior, ::AdvancedHMC.AbstractHMCSampler)
AbstractMCMC.sample(rng::AbstractRNG, ::Comrade.VLBIPosterior, ::NestedSamplers.Nested)
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
```

## Internal (Not Public API)

### eht-imaging interface (Internal)

```@docs
Comrade.extract_amp
Comrade.extract_cphase
Comrade.extract_lcamp
Comrade.extract_vis
Comrade.extract_coherency
```


