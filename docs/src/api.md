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
Comrade.caltable(::Comrade.SiteLookup, array::AbstractArray)
Comrade.caltable(::Comrade.SiteArray)
Comrade.JonesSandwich
Comrade.VLBIModel
Comrade.CalPrior
Comrade.CalPrior(::NamedTuple, ::JonesCache)
Comrade.CalPrior(::NamedTuple, ::NamedTuple, ::JonesCache)
Comrade.RIMEModel
Comrade.ObsSegmentation
Comrade.IntegSeg
Comrade.ScanSeg
Comrade.TrackSeg
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
Comrade.domain
Comrade.timetable
Comrade.sites
Comrade.uvpositions
Comrade.ArrayConfiguration
Comrade.ClosureConfig
Comrade.AbstractInterferometryDatum
Comrade.EHTArrayBaselineDatum
Comrade.EHTObservationTable
Comrade.EHTArrayConfiguration
Comrade.EHTCoherencyDatum
Comrade.EHTVisibilityDatum
Comrade.EHTVisibilityAmplitudeDatum
Comrade.EHTLogClosureAmplitudeDatum
Comrade.EHTClosurePhaseDatum
Comrade.Scan
Comrade.TimeTable
```

## Fourier Transform Algorithms
```@docs
VLBISkyModels.NFFTAlg()
VLBISkyModels.DFTAlg()
VLBISkyModels.FFTAlg()
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


