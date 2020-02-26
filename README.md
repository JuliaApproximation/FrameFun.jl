
[![Build Status](https://travis-ci.org/JuliaApproximation/FrameFun.jl.svg?branch=master)](https://travis-ci.org/JuliaApproximation/BasisFunctions.jl.svg?branch=master)
[![Coverage Status](https://coveralls.io/repos/github/JuliaApproximation/FrameFun.jl/badge.svg)](https://coveralls.io/github/JuliaApproximation/FrameFun.jl)

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
Plots.savefig(P,"images/lowprecision.png")
```

![](images/lowprecision.png)

The bases support any AbstractFloat subtype, so high precision approximations are straightforward:



```julia
B = Fourier(61) → big(-1)..big(1)
F = Fun(f,B,D)

P = plot(F,layout=2)
plot!(F,f,subplot=2)
Plots.savefig(P,"images/highprecision.png")
```

![](images/highprecision.png)

![](images/highprecision.png)

# Frame Approximations in 2D

In higher dimensions, a basis can be any tensorproduct of (scaled) lower dimensional bases:


```julia
using StaticArrays
C = disk(1.0)\disk(0.3,SVector(0.2, 0.5))
B = (Fourier(31) → -1.3..1.3)^2
f = (x,y)->exp(x+y)
F = Fun(f,B,C)

P = heatmap(F,layout=2,aspect_ratio=1)
plot!(F,f,subplot=2,aspect_ratio=1)
Plots.savefig(P,"images/deathstar.png")
```

![](images/deathstar.png)

Even fractal domains are not a problem:


```julia
B = (Fourier(31) → -1.0..0.35) ⊗ (Fourier(31) → -0.65..0.65)
f = (x,y)->cos(10*x*y)
F = Fun(f, B, mandelbrot())

P = heatmap(F,layout=2,aspect_ratio=1)
plot!(F,f,aspect_ratio=1,subplot=2)
Plots.savefig(P,"images/mandelbrot")
```

    GKS: Rectangle definition is invalid in routine SET_WINDOW
    GKS: Rectangle definition is invalid in routine CELLARRAY
    origin outside current window


![](images/mandelbrot.png)


## Installation

FrameFun.jl is not added to the Julia General registry and depends on the unregistered packages GridArrays.jl and BasisFunctions.jl

### Recommended
For Julia 1.1 or higher, you can add the FrameFun registry.
From the Julia REPL, type `]` to enter Pkg mode and run

```julia
pkg> registry add https://github.com/FramefunVC/FrameFunRegistry
pkg> add FrameFun
```

### Legacy
In Julia 1.0, the packages can be installed by cloning their git repository. From the Julia REPL, type `]` to enter Pkg mode and run

```julia
pkg> add https://github.com/JuliaApproximation/BasisFunctions.jl
pkg> add https://github.com/JuliaApproximation/FrameFun.jl
```
