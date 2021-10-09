
| **Build Status** | **Coverage** |
|------------------|--------------|
| [![Build Status](https://github.com/JuliaApproximation/FrameFun.jl/workflows/CI/badge.svg?branch=master)](https://github.com/JuliaApproximation/FrameFun.jl/actions/workflows/CI.yml) | [![Coverage](https://codecov.io/gh/JuliaApproximation/FrameFun.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaApproximation/FrameFun.jl)

FrameFun
========

Exploring practical possibilities of approximating functions with frames rather than with a basis. The package is heavily inspired by the Chebfun project and the Julia package ApproxFun.


```julia
using BasisFunctions, Plots, DomainSets, FrameFun
gr();
```

# Frame Approximations in 1D

After choosing a suitable Basis and Domain, any function can be approximated in the resulting frame:


```julia
B = Fourier(61) → -1..1
D = -0.5..0.5
f = x->x
F = Fun(f,B,D)

P = plot(F,layout = 2)
plot!(F,f,subplot=2)
```

![](images/lowprecision.png)

The bases support any AbstractFloat subtype, so high precision approximations are straightforward:



```julia
B = Fourier(61) → big(-1)..big(1)
F = Fun(f,B,D)

P = plot(F,layout=2)
plot!(F,f,subplot=2)
```

![](images/highprecision.png)

![](images/highprecision.png)

# Frame Approximations in 2D

In higher dimensions, a basis can be any tensorproduct of (scaled) lower dimensional bases:


```julia
using StaticArrays
C = Disk(1.0)\Disk(0.3,SVector(0.2, 0.5))
B = (Fourier(31) → -1.3..1.3)^2
f = (x,y)->exp(x+y)
F = Fun(f,B,C)

P = heatmap(F,layout=2,aspect_ratio=1)
plot!(F,f,subplot=2,aspect_ratio=1)
```

![](images/deathstar.png)

Even fractal domains are not a problem:


```julia
B = (Fourier(31) → -1.0..0.35) ⊗ (Fourier(31) → -0.65..0.65)
f = (x,y)->cos(10*x*y)
F = Fun(f, B, mandelbrot())

P = heatmap(F,layout=2,aspect_ratio=1)
plot!(F,f,aspect_ratio=1,subplot=2)
```

![](images/mandelbrot.png)


## Installation

From the Julia REPL, type `]` to enter Pkg mode and run

```julia
pkg> add DomainSets, BasisFunctions, FrameFun
```
