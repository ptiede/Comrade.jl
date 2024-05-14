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

### Calibration Models

```@docs
Comrade.CalTable
Comrade.caltable(::Comrade.SiteArray)
Comrade.IIDSitePrior
Comrade.ArrayPrior
Comrade.RIMEModel
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
Comrade.IdealInstrument
Comrade.InstrumentModel
Comrade.site_tuple
Comrade.SiteArray
Comrade.SiteLookup
```

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



## Bayesian Tools

### Posterior Constructions

```@docs
Comrade.ascube
Comrade.asflat
Comrade.flatten
Comrade.inverse
Comrade.prior_sample
Comrade.likelihood
Comrade.simulate_observation
Comrade.dataproducts
Comrade.skymodel
Comrade.instrumentmodel
Comrade.vlbimodel
Comrade.sample(::Posterior)
Comrade.transform
Comrade.MultiRadioLikelihood
Comrade.Posterior
Comrade.TransformedPosterior
Comrade.RadioLikelihood
Comrade.IsFlat
Comrade.IsCube
```

### Sampler Tools
```@docs
Comrade.samplertype
```

## Misc

```@docs
Comrade.site_tuple
Comrade.dirty_image
Comrade.dirty_beam
Comrade.beamsize
```

## Internal (Not Public API)
```@docs
Comrade.extract_FRs
ComradeBase._visibilitymap!
ComradeBase._visibilitymap
```

### eht-imaging interface (Internal)

```@docs
Comrade.extract_amp
Comrade.extract_cphase
Comrade.extract_lcamp
Comrade.extract_vis
Comrade.extract_coherency
```


