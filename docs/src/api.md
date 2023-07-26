# Comrade API



## Contents

```@contents
Pages = ["api.md"]
```

## Index

```@index
Pages = ["api.md"]
```

## Model Definitions

### Calibration Models

```@docs
Comrade.corrupt
Comrade.CalTable
Comrade.caltable(::Comrade.JonesCache, ::AbstractVector)
Comrade.caltable(::Comrade.EHTObservation, ::AbstractVector)
Comrade.DesignMatrix
Comrade.JonesCache
Comrade.TransformCache
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
Comrade.jonesT
Base.map(::Any, ::Vararg{Comrade.JonesPairs})
Comrade.PoincareSphere2Map
Comrade.caltable
Comrade.JonesPairs
Comrade.GainSchema
```

### Models

```@docs
VLBISkyModels.DFTAlg(::Comrade.ArrayConfiguration)
VLBISkyModels.DFTAlg(::Comrade.EHTObservation)
VLBISkyModels.NFFTAlg(::Comrade.ArrayConfiguration)
VLBISkyModels.NFFTAlg(::Comrade.EHTObservation)
```


### Polarized Models

```@docs
PolarizedTypes.mbreve
PolarizedTypes.m̆
PolarizedTypes.evpa
PolarizedTypes.linearpol
```



## Data Types


```@docs
Comrade.extract_table
Comrade.ComplexVisibilities
Comrade.VisibilityAmplitudes
Comrade.ClosurePhases
Comrade.LogClosureAmplitudes
Comrade.Coherencies
Comrade.amplitude(::Comrade.EHTVisibilityDatum)
Comrade.amplitude(::Comrade.EHTVisibilityAmplitudeDatum)
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
Comrade.EHTClosurePhaseDatum
Comrade.EHTLogClosureAmplitudeDatum
Comrade.EHTVisibilityDatum
Comrade.EHTVisibilityAmplitudeDatum
Comrade.Scan
Comrade.ScanTable
```

## eht-imaging interface (Internal)

```@docs
Comrade.extract_amp
Comrade.extract_cphase
Comrade.extract_lcamp
Comrade.extract_vis
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

## Misc

```@docs
Comrade.station_tuple
Comrade.dirty_image
Comrade.dirty_beam
```

## Internal (Not Public API)
```@docs
Comrade.extract_FRs
```

