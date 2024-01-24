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
Comrade.caltable(::Comrade.JonesCache, ::AbstractVector)
Comrade.caltable(::Comrade.EHTObservation, ::AbstractVector)
Comrade.DesignMatrix
Comrade.JonesCache
Comrade.ResponseCache
Comrade.JonesModel
Comrade.VLBIModel
Comrade.CalPrior
Comrade.CalPrior(::NamedTuple, ::JonesCache)
Comrade.CalPrior(::NamedTuple, ::NamedTuple, ::JonesCache)
Comrade.RIMEModel
Comrade.ObsSegmentation
Comrade.IntegSeg
Comrade.ScanSeg
Comrade.TrackSeg
Comrade.FixedSeg
Comrade.jonescache(::Comrade.EHTObservation, ::Comrade.ObsSegmentation)
Comrade.SingleReference
Comrade.RandomReference
Comrade.SEFDReference
Comrade.jonesStokes
Comrade.jonesG
Comrade.jonesD
Comrade.jonesR
Base.map(::Any, ::Vararg{Comrade.JonesPairs})
Comrade.caltable
Comrade.JonesPairs
Comrade.GainSchema
Comrade.SegmentedJonesCache
```

### Models

For the description of the model API see [VLBISkyModels](https://ehtjulia.github.io/VLBISkyModels.jl/stable/).




## Data Types


```@docs
Comrade.extract_table
Comrade.ComplexVisibilities
Comrade.VisibilityAmplitudes
Comrade.ClosurePhases
Comrade.LogClosureAmplitudes
Comrade.Coherencies
Comrade.baselines
Comrade.arrayconfig
Comrade.closure_phase(::Comrade.EHTVisibilityDatum, ::Comrade.EHTVisibilityDatum, ::Comrade.EHTVisibilityDatum)
Comrade.getdata
Comrade.getuv
Comrade.getuvtimefreq
Comrade.scantable
Comrade.stations
Comrade.uvpositions
Comrade.ArrayConfiguration
Comrade.ClosureConfig
Comrade.AbstractInterferometryDatum
Comrade.ArrayBaselineDatum
Comrade.EHTObservation
Comrade.EHTArrayConfiguration
Comrade.EHTCoherencyDatum
Comrade.EHTVisibilityDatum
Comrade.EHTVisibilityAmplitudeDatum
Comrade.EHTLogClosureAmplitudeDatum
Comrade.EHTClosurePhaseDatum
Comrade.Scan
Comrade.ScanTable
```

## Model Cache
```@docs
VLBISkyModels.NFFTAlg(::Comrade.EHTObservation)
VLBISkyModels.NFFTAlg(::Comrade.ArrayConfiguration)
VLBISkyModels.DFTAlg(::Comrade.EHTObservation)
VLBISkyModels.DFTAlg(::Comrade.ArrayConfiguration)
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
Comrade.station_tuple
Comrade.dirty_image
Comrade.dirty_beam
Comrade.beamsize
```

## Internal (Not Public API)
```@docs
Comrade.extract_FRs
ComradeBase._visibilities!
ComradeBase._visibilities
```

### eht-imaging interface (Internal)

```@docs
Comrade.extract_amp
Comrade.extract_cphase
Comrade.extract_lcamp
Comrade.extract_vis
Comrade.extract_coherency
```


