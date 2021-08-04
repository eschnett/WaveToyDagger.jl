# WaveToyDagger.jl

A WaveToy based on [Dagger.jl](https://github.com/JuliaParallel/Dagger.jl).

```julia
julia> using Revise

julia> using Distributed

julia> addprocs(6, exeflags="--project");

julia> @everywhere using WaveToyDagger
[ Info: Precompiling WaveToyDagger [1ebf5abb-aa52-4892-b9b8-f6b456a24862]

julia> using BenchmarkTools

julia> @benchmark WaveToyDagger.main()
BenchmarkTools.Trial: 25 samples with 1 evaluation.
 Range (min … max):  156.787 ms … 334.186 ms  ┊ GC (min … max): 0.00% … 14.22%
 Time  (median):     182.191 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   208.185 ms ±  53.668 ms  ┊ GC (mean ± σ):  4.29% ±  6.53%

   █   ▁
  ▆█▆▆▆█▆▆▆▆▁▁▁▆▁▁▆▁▁▁▆▁▁▁▁▁▆▁▁▆▆▁▁▁▁▁▁▁▁▁▆▁▁▆▁▁▁▁▆▁▁▆▁▁▁▁▁▁▁▁▆ ▁
  157 ms           Histogram: frequency by time          334 ms <

 Memory estimate: 10.38 MiB, allocs estimate: 198917.
```
