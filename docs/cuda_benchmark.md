# Project Benchmark Results (CPU vs GPU)

This file records the benchmark results for both CPU and GPU. Please note that the CPU related tests are conducted on `Float64` type while the GPU related tests are conducted on `Float32` type.

1. Benchmarks on `elixir_advection_basic.jl`
```Julia
# Benchmark result for `advection_basic_1d.jl` on CPU

BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  223.597 μs …  4.951 ms  ┊ GC (min … max): 0.00% … 94.27%
 Time  (median):     228.619 μs              ┊ GC (median):    0.00%
 Time  (mean ± σ):   231.565 μs ± 66.704 μs  ┊ GC (mean ± σ):  0.40% ±  1.33%

      ▁▅██▅▁                                                    
  ▂▃▅▆██████▇▆▆▆▅▅▄▄▄▄▃▃▃▂▂▂▂▂▂▂▃▃▃▃▃▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▁▂▂▂▂ ▃
  224 μs          Histogram: frequency by time          260 μs <

 Memory estimate: 12.17 KiB, allocs estimate: 51.

# Benchmark result for `advection_basic_1d.jl` on GPU

BenchmarkTools.Trial: 47 samples with 1 evaluation.
 Range (min … max):   87.600 ms … 219.162 ms  ┊ GC (min … max): 0.00% … 6.80%
 Time  (median):      90.206 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   106.622 ms ±  40.771 ms  ┊ GC (mean ± σ):  1.98% ± 2.65%

  ▄█                                                             
  ██▃▁▁▅▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▅▃▁▁▃ ▁
  87.6 ms          Histogram: frequency by time          219 ms <

 Memory estimate: 6.33 MiB, allocs estimate: 112029.
```

```Julia
# Benchmark result for `advection_basic_2d.jl` on CPU

BenchmarkTools.Trial: 1175 samples with 1 evaluation.
 Range (min … max):  4.221 ms …   9.192 ms  ┊ GC (min … max): 0.00% … 52.16%
 Time  (median):     4.240 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   4.249 ms ± 145.701 μs  ┊ GC (mean ± σ):  0.10% ±  1.52%

        ▁▃▄▆▆█▄▆▄▃                                             
  ▂▃▃▅▆▇██████████▇▇▆█▆▆▅▅▄▄▄▄▃▃▃▃▃▃▃▄▄▃▃▃▂▃▃▄▂▃▂▂▂▁▂▁▂▂▁▂▂▂▂ ▄
  4.22 ms         Histogram: frequency by time        4.31 ms <

 Memory estimate: 72.00 KiB, allocs estimate: 53.
 
# Benchmark result for `advection_basic_2d.jl` on GPU

BenchmarkTools.Trial: 2363 samples with 1 evaluation.
 Range (min … max):  1.626 ms … 46.685 ms  ┊ GC (min … max): 0.00% … 20.89%
 Time  (median):     1.756 ms              ┊ GC (median):    0.00%
 Time  (mean ± σ):   2.103 ms ±  3.831 ms  ┊ GC (mean ± σ):  3.40% ±  1.78%

         ▁▁▄▂▁▁▃▇▇█▅▅▇▄▄▃▁▁ ▁                                 
  ▂▂▃▃▃▃▅████████████████████▆▅▅▆▄▄▃▃▄▃▃▃▂▄▂▂▂▂▂▃▂▂▂▂▁▂▂▁▁▁▁ ▄
  1.63 ms        Histogram: frequency by time        2.03 ms <

 Memory estimate: 345.16 KiB, allocs estimate: 2056.
```

```Julia
# Benchmark result for `advection_basic_3d.jl` on CPU

BenchmarkTools.Trial: 17 samples with 1 evaluation.
 Range (min … max):  301.413 ms … 303.714 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     301.967 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   302.104 ms ± 573.530 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

               ▃▃   █                                            
  ▇▁▁▁▇▁▁▇▇▁▁▁▇██▇▇▁█▁▁▁▁▁▇▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▇▁▁▁▁▁▁▁▁▁▁▇ ▁
  301 ms           Histogram: frequency by time          304 ms <

 Memory estimate: 3.26 MiB, allocs estimate: 66.

# Benchmark result for `advection_basic_3d.jl` on GPU

BenchmarkTools.Trial: 12 samples with 1 evaluation.
 Range (min … max):  433.459 ms … 441.928 ms  ┊ GC (min … max): 13.32% … 13.48%
 Time  (median):     435.328 ms               ┊ GC (median):    13.32%
 Time  (mean ± σ):   435.897 ms ±   2.281 ms  ┊ GC (mean ± σ):  13.35% ±  0.10%

  ▁   ▁▁   ▁ ▁ █  ▁  ▁▁              ▁                        ▁  
  █▁▁▁██▁▁▁█▁█▁█▁▁█▁▁██▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  433 ms           Histogram: frequency by time          442 ms <

 Memory estimate: 839.19 MiB, allocs estimate: 151753.
```

2. Benchmarks on `elixir_euler_ec.jl`
```Julia
# Benchmark result for `euler_ec_1d.jl` on CPU

BenchmarkTools.Trial: 498 samples with 1 evaluation.
 Range (min … max):   9.791 ms …  10.793 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     10.114 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   10.031 ms ± 144.959 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

         ▂▂▁                                   ▁▂▇█▄▁▂          
  ▃▃▄▄▆▇▇███▅▄▆▅▃▃▃▃▃▂▃▂▂▂▁▂▂▁▁▁▁▁▁▁▁▁▁▁▁▃▂▃▆▇████████▅▆▃▃▄▄▃▃ ▃
  9.79 ms         Histogram: frequency by time         10.2 ms <

 Memory estimate: 55.08 KiB, allocs estimate: 53.
 
# Benchmark result for `euler_ec_1d.jl` on GPU

BenchmarkTools.Trial: 13 samples with 1 evaluation.
 Range (min … max):  347.788 ms … 422.485 ms  ┊ GC (min … max): 0.00% … 2.51%
 Time  (median):     419.465 ms               ┊ GC (median):    2.72%
 Time  (mean ± σ):   410.706 ms ±  23.081 ms  ┊ GC (mean ± σ):  2.35% ± 1.02%

                                                            █ ▂  
  ▅▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▅▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▅█▁█ ▁
  348 ms           Histogram: frequency by time          422 ms <

 Memory estimate: 39.82 MiB, allocs estimate: 407391.
```

```Julia
# Benchmark result for `euler_ec_2d.jl` on CPU

BenchmarkTools.Trial: 3 samples with 1 evaluation.
 Range (min … max):  2.495 s …  2.497 s  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     2.497 s             ┊ GC (median):    0.00%
 Time  (mean ± σ):   2.496 s ± 1.366 ms  ┊ GC (mean ± σ):  0.00% ± 0.00%

  █                                              █       █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁█ ▁
  2.49 s        Histogram: frequency by time         2.5 s <

 Memory estimate: 8.01 MiB, allocs estimate: 69.
 
# Benchmark result for `euler_ec_2d.jl` on GPU

BenchmarkTools.Trial: 3 samples with 1 evaluation.
 Range (min … max):  2.117 s …  2.124 s  ┊ GC (min … max): 17.46% … 17.40%
 Time  (median):     2.120 s             ┊ GC (median):    17.44%
 Time  (mean ± σ):   2.120 s ± 3.398 ms  ┊ GC (mean ± σ):  17.44% ±  0.04%

  █                    █                                 █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  2.12 s        Histogram: frequency by time        2.12 s <

 Memory estimate: 4.12 GiB, allocs estimate: 483170.
```

```Julia
# Benchmark result for `euler_ec_3d.jl` on CPU

BenchmarkTools.Trial: 3 samples with 1 evaluation.
 Range (min … max):  2.394 s …  2.400 s  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     2.399 s             ┊ GC (median):    0.00%
 Time  (mean ± σ):   2.398 s ± 3.421 ms  ┊ GC (mean ± σ):  0.00% ± 0.00%

  █                                              █       █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁█ ▁
  2.39 s        Histogram: frequency by time         2.4 s <

 Memory estimate: 20.01 MiB, allocs estimate: 68.

# Benchmark result for `euler_ec_3d.jl` on GPU

BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.447 s …  1.455 s  ┊ GC (min … max): 12.75% … 12.76%
 Time  (median):     1.450 s             ┊ GC (median):    12.76%
 Time  (mean ± σ):   1.450 s ± 3.366 ms  ┊ GC (mean ± σ):  12.76% ±  0.01%

  █   █                         █                        █  
  █▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.45 s        Histogram: frequency by time        1.45 s <

 Memory estimate: 3.69 GiB, allocs estimate: 136804.
```

3. Benchmarks on `elixir_euler_source_terms.jl`
```Julia
# Benchmark result for `euler_source_terms_1d.jl` on CPU

BenchmarkTools.Trial: 1514 samples with 1 evaluation.
 Range (min … max):  3.252 ms …   7.680 ms  ┊ GC (min … max): 0.00% … 56.65%
 Time  (median):     3.292 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   3.294 ms ± 114.766 μs  ┊ GC (mean ± σ):  0.09% ±  1.46%

             ▂▂▂▁    ▃▄█▅▆▅▃                                   
  ▂▂▃▃▃▄▃▃▅▆█████▆▆██████████▇▄▄▄▃▃▃▃▃▂▃▂▂▃▂▂▂▂▁▁▂▁▁▁▂▂▂▂▂▂▁▂ ▄
  3.25 ms         Histogram: frequency by time        3.36 ms <

 Memory estimate: 37.03 KiB, allocs estimate: 54.
 
# Benchmark result for `euler_source_terms_1d.jl` on GPU

BenchmarkTools.Trial: 17 samples with 1 evaluation.
 Range (min … max):  245.944 ms … 369.302 ms  ┊ GC (min … max): 0.00% … 3.55%
 Time  (median):     338.529 ms               ┊ GC (median):    3.87%
 Time  (mean ± σ):   300.392 ms ±  48.204 ms  ┊ GC (mean ± σ):  2.52% ± 2.15%

  ▄ ▁                                          ▁█                
  █▆█▆▁▁▁▁▁▁▁▁▆▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁██▆▆▁▁▁▁▁▁▁▁▁▁▁▆ ▁
  246 ms           Histogram: frequency by time          369 ms <

 Memory estimate: 23.15 MiB, allocs estimate: 319021.
```

```Julia
# Benchmark result for `euler_source_terms_2d.jl` on CPU

BenchmarkTools.Trial: 16 samples with 1 evaluation.
 Range (min … max):  311.442 ms … 316.894 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     312.043 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   312.675 ms ±   1.427 ms  ┊ GC (mean ± σ):  0.00% ± 0.00%

  ▃   █                                                          
  █▁▁▇█▇▇▇▁▁▇▇▁▁▁▁▁▇▁▁▁▁▁▇▁▁▇▁▁▁▁▁▁▇▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▇ ▁
  311 ms           Histogram: frequency by time          317 ms <

 Memory estimate: 2.01 MiB, allocs estimate: 70.
 
# Benchmark result for `euler_source_terms_2d.jl` on GPU

BenchmarkTools.Trial: 7 samples with 1 evaluation.
 Range (min … max):  758.963 ms … 771.515 ms  ┊ GC (min … max): 11.07% … 11.68%
 Time  (median):     765.217 ms               ┊ GC (median):    11.30%
 Time  (mean ± σ):   765.419 ms ±   4.784 ms  ┊ GC (mean ± σ):  11.34% ±  0.21%

  █     █                   █   █      █                    █ █  
  █▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁█▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁█ ▁
  759 ms           Histogram: frequency by time          772 ms <

 Memory estimate: 899.09 MiB, allocs estimate: 440106.
```

```Julia
# Benchmark result for `euler_source_terms_3d.jl` on CPU

BenchmarkTools.Trial: 11 samples with 1 evaluation.
 Range (min … max):  488.843 ms … 490.704 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     489.669 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   489.703 ms ± 528.316 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

  ▁              ▁ ▁   █     ▁ ▁█                        ▁    ▁  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁█▁▁▁█▁▁▁▁▁█▁██▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁█ ▁
  489 ms           Histogram: frequency by time          491 ms <

 Memory estimate: 2.51 MiB, allocs estimate: 68.
 
# Benchmark result for `euler_source_terms_3d.jl` on GPU

BenchmarkTools.Trial: 7 samples with 1 evaluation.
 Range (min … max):  773.680 ms … 782.401 ms  ┊ GC (min … max): 8.68% … 9.31%
 Time  (median):     779.970 ms               ┊ GC (median):    9.36%
 Time  (mean ± σ):   779.598 ms ±   2.801 ms  ┊ GC (mean ± σ):  9.26% ± 0.26%

  ▁                                     ▁    █      ▁▁        ▁  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁█▁▁▁▁▁▁██▁▁▁▁▁▁▁▁█ ▁
  774 ms           Histogram: frequency by time          782 ms <

 Memory estimate: 1.21 GiB, allocs estimate: 387978.
```
