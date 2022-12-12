# Benchmarks

## Imaging no gains

### Comrade 0.7.0

```
julia> @benchmark ℓ($(rand(ndim)))
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):   9.870 μs …  14.060 ms  ┊ GC (min … max): 0.00% … 96.30%
 Time  (median):     12.170 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   14.064 μs ± 140.647 μs  ┊ GC (mean ± σ):  9.63% ±  0.96%

       ▁▆███▇▇▇▆▃                                               
  ▁▁▁▃▅███████████▇▅▃▂▂▂▁▁▁▁▁▁▁▁▁▁▂▂▂▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▃
  9.87 μs         Histogram: frequency by time         22.6 μs <

 Memory estimate: 16.67 KiB, allocs estimate: 11.
```

```
julia> @benchmark Zygote.gradient($ℓ, $(rand(ndim)))
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  141.594 μs …  10.580 ms  ┊ GC (min … max):  0.00% … 95.36%
 Time  (median):     151.444 μs               ┊ GC (median):     0.00%
 Time  (mean ± σ):   212.204 μs ± 499.147 μs  ┊ GC (mean ± σ):  15.32% ±  6.63%

     ▅█                                                          
  ▁▂▆███▄▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▂▃▄▄▄▃▃▂▂▁▁▁ ▂
  142 μs           Histogram: frequency by time          263 μs <

 Memory estimate: 917.08 KiB, allocs estimate: 966.

julia> 
```

### Comrade 0.6.9

```
julia> @benchmark ℓ($(rand(ndim)))
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  14.750 μs …  20.074 ms  ┊ GC (min … max):  0.00% … 99.90%
 Time  (median):     17.450 μs               ┊ GC (median):     0.00%
 Time  (mean ± σ):   20.011 μs ± 200.726 μs  ┊ GC (mean ± σ):  10.02% ±  1.00%

    ▆█▁  ▂▄▆▆▅▆▇▅▂                                              
  ▂▇███▅▇█████████▆▅▄▃▂▂▂▂▂▂▁▂▂▂▂▂▂▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▃
  14.8 μs         Histogram: frequency by time         29.9 μs <

 Memory estimate: 25.83 KiB, allocs estimate: 16.
```

```
julia> @benchmark Zygote.gradient($ℓ, $(rand(ndim)))
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  250.446 μs …  10.622 ms  ┊ GC (min … max):  0.00% … 93.77%
 Time  (median):     270.746 μs               ┊ GC (median):     0.00%
 Time  (mean ± σ):   321.700 μs ± 514.504 μs  ┊ GC (mean ± σ):  11.67% ±  7.01%

    ▂▄▅▅▇██▇▇▆▅▄▃▃▂▂▁                         ▁▁▂▂▂▃▂▂▂▂▂▂▂▁▁   ▃
  ▅▇████████████████████▇▆▇▅▅▅▅▆▄▄▅▅▅▁▃▁▄▄▃▆▆█████████████████▆ █
  250 μs        Histogram: log(frequency) by time        391 μs <

 Memory estimate: 1.06 MiB, allocs estimate: 3612.
```

## Imaging with gains 

Fitting Amplitudes and Closure Phases (diagonal)

### Comrade 0.7.0-dev

```
julia> @benchmark ℓ($(rand(ndim)))
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  12.660 μs …  21.255 ms  ┊ GC (min … max):  0.00% … 97.91%
 Time  (median):     15.570 μs               ┊ GC (median):     0.00%
 Time  (mean ± σ):   17.741 μs ± 212.506 μs  ┊ GC (mean ± σ):  11.73% ±  0.98%

       ▃▃▃▁        ▃▅▅██▇▅▄▃▁                                   
  ▁▂▃▅█████▇▅▄▃▃▄▆█████████████▆▅▄▄▃▃▃▂▂▂▂▂▂▂▂▁▂▂▁▁▁▁▁▁▁▁▁▁▂▁▂ ▄
  12.7 μs         Histogram: frequency by time         21.1 μs <

 Memory estimate: 29.98 KiB, allocs estimate: 26.
```

```
julia> @benchmark Zygote.gradient($ℓ, $(rand(ndim)))
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  195.526 μs …  12.042 ms  ┊ GC (min … max):  0.00% … 93.48%
 Time  (median):     209.486 μs               ┊ GC (median):     0.00%
 Time  (mean ± σ):   277.949 μs ± 591.254 μs  ┊ GC (mean ± σ):  16.42% ±  7.67%

   ▃▇██▆▄▃▂▂▁▁                                ▂▂▁ ▁▂▃▂▂▂▂▁▁     ▂
  ▇█████████████▆▆▆▃▄▅▅▅▃▄▁▄▄▁▃▃▄▃▃▄▃▃▁▁▃▁▄▄▄▇███████████████▆▅ █
  196 μs        Histogram: log(frequency) by time        398 μs <

 Memory estimate: 1.50 MiB, allocs estimate: 1180.
```


### Comrade 0.6.8

```
julia> @benchmark ℓ($(rand(ndim)))
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  16.469 μs …  18.204 ms  ┊ GC (min … max):  0.00% … 97.49%
 Time  (median):     19.919 μs               ┊ GC (median):     0.00%
 Time  (mean ± σ):   23.538 μs ± 256.765 μs  ┊ GC (mean ± σ):  15.08% ±  1.38%

     ▅▇█▇▃        ▂▄▇▆▆▆▄▃▂                                     
  ▁▂▅██████▄▂▂▂▃▄▇██████████▇▆▄▄▃▃▃▂▃▃▂▂▂▂▂▂▁▁▁▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▃
  16.5 μs         Histogram: frequency by time         27.6 μs <

 Memory estimate: 37.98 KiB, allocs estimate: 23.
```

```
julia> @benchmark Zygote.gradient($ℓ, $(rand(ndim)))
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  307.953 μs …   6.769 ms  ┊ GC (min … max):  0.00% … 89.39%
 Time  (median):     327.218 μs               ┊ GC (median):     0.00%
 Time  (mean ± σ):   389.707 μs ± 569.119 μs  ┊ GC (mean ± σ):  13.56% ±  8.72%

  ▂▆███▇▆▅▄▃▃▂▁▁▁▂▁                                             ▂
  ██████████████████▇▇▆▆▆▆▆▅▅▆▄▅▆▆▅▆▆▇▇▆▆▇▇▇▆▆▇▆▇▆▆▄▅▄▄▄▄▄▁▃▁▁▄ █
  308 μs        Histogram: log(frequency) by time        583 μs <

 Memory estimate: 1.64 MiB, allocs estimate: 3805.
```

## M-Ring fitting (paper_example.jl)


### Comrade 0.7.0
```
julia> @benchmark ℓ($(rand(ndim)))
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  30.819 μs …   8.158 ms  ┊ GC (min … max): 0.00% … 90.67%
 Time  (median):     33.910 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   37.833 μs ± 156.955 μs  ┊ GC (mean ± σ):  7.70% ±  1.85%

         ▁▂▆█▆▄▂▁                                               
  ▁▂▃▅▆██████████▇▆▆▇▇▆▆▆▆▅▄▅▄▄▃▃▂▂▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▃
  30.8 μs         Histogram: frequency by time         45.1 μs <

 Memory estimate: 91.48 KiB, allocs estimate: 49.
```

```
julia> @benchmark gf($(rand(10)))
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  191.656 μs …   7.624 ms  ┊ GC (min … max):  0.00% … 95.59%
 Time  (median):     206.831 μs               ┊ GC (median):     0.00%
 Time  (mean ± σ):   284.461 μs ± 572.746 μs  ┊ GC (mean ± σ):  18.09% ±  8.71%

  ▄▅██▅▄▃▂                             ▁▂▂▂▁▁  ▂▂▁      ▁▁▂▁    ▂
  █████████▇▆▄▁▄▃▃▃▃▁▄▃▁▁▁▁▁▁▁▁▁▁▄▆▆▅▄▇██████▆▇█████▆▅▆▇█████▇▇ █
  192 μs        Histogram: log(frequency) by time        486 μs <

 Memory estimate: 1.92 MiB, allocs estimate: 351.
```


### Comrade 0.6.8
```
julia> @benchmark ℓ($(rand(ndim)))
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  25.480 μs … 566.369 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     27.430 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   27.712 μs ±   5.652 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

              ▁▇██▆▄▄▃▁                                         
  ▁▁▁▁▁▁▁▂▄██▇██████████▅▄▃▃▂▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▃
  25.5 μs         Histogram: frequency by time         32.7 μs <

 Memory estimate: 12.45 KiB, allocs estimate: 10.
```


```
julia> gf = Comrade.make_pullback(ℓ, AD.ForwardDiffBackend{10}())
julia> @benchmark gf($(rand(ndim)))
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  74.049 μs …   6.981 ms  ┊ GC (min … max): 0.00% … 96.26%
 Time  (median):     75.868 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   80.779 μs ± 151.681 μs  ┊ GC (mean ± σ):  4.09% ±  2.15%

   ▂▅▇██▇▆▅▃▃▂▂▁  ▁▁▁                         ▁▁▁▁▁▁▁   ▁▁     ▂
  ▇████████████████████▇▇▆▅▅▅▆▄▅▅▆▄▄▃▅▄▅▅▄▆▇████████████████▇▇ █
  74 μs         Histogram: log(frequency) by time      93.3 μs <

 Memory estimate: 130.22 KiB, allocs estimate: 17.
```


# Computing Environment

```
Julia Version 1.8.3
Commit 0434deb161e (2022-11-14 20:14 UTC)
Platform Info:
  OS: Linux (x86_64-linux-gnu)
  CPU: 32 × AMD Ryzen 9 7950X 16-Core Processor
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-13.0.1 (ORCJIT, znver3)
  Threads: 1 on 32 virtual cores
Environment:
  JULIA_EDITOR = code
  JULIA_NUM_THREADS = 1
```