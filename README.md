
[![Build Status](https://travis-ci.org/daanhb/FrameFun.jl.svg?branch=master)](https://travis-ci.org/daanhb/FrameFun.jl)
[![Coverage Status](https://coveralls.io/repos/github/daanhb/FrameFun.jl/badge.svg)](https://coveralls.io/github/daanhb/FrameFun.jl)

FrameFun
========

Exploring practical possibilities of approximating functions with frames rather than with a basis. The package is heavily inspired by the Chebfun project and the Julia package ApproxFun.

# Frame Approximations in 1D

After choosing a suitable Basis and Domain, any function can be approximated in the resulting frame:


```julia
using BasisFunctions
using Plots;gr()
using Domains
using FrameFun
B = FourierBasis(61, -1, 1)
D = interval(-0.5,0.5)
f = x->x
F = Fun(f,B,D)

P = plot(F,plot_ext=true, layout = 2)
plot!(F,f,plot_ext=true, subplot=2)
Plots.savefig(P,"images/lowprecision.png")
```

    WARNING: Method definition infimum(Domains.Domain{T} where T) in module Domains at /Users/vincentcp/julia/Domains/src/generic/domain.jl:84 overwritten in module FrameFun at /Users/vincentcp/julia/FrameFun/src/domains/extensions.jl:86.
    WARNING: Method definition supremum(Domains.Domain{T} where T) in module Domains at /Users/vincentcp/julia/Domains/src/generic/domain.jl:85 overwritten in module FrameFun at /Users/vincentcp/julia/FrameFun/src/domains/extensions.jl:88.
    WARNING: Method definition infimum(Domains.ProductDomain{DD, S, T} where T where S where DD) in module Domains at /Users/vincentcp/julia/Domains/src/generic/productdomain.jl:108 overwritten in module FrameFun at /Users/vincentcp/julia/FrameFun/src/domains/extensions.jl:90.
    WARNING: Method definition supremum(Domains.ProductDomain{DD, S, T} where T where S where DD) in module Domains at /Users/vincentcp/julia/Domains/src/generic/productdomain.jl:109 overwritten in module FrameFun at /Users/vincentcp/julia/FrameFun/src/domains/extensions.jl:92.


![](images/lowprecision.png)

The bases support any AbstractFloat subtype, so high precision approximations are straightforward:



```julia
B = FourierBasis(61, BigFloat(-1), BigFloat(1))
F = Fun(f,B,D)

P = plot(F,plot_ext=true,layout=2)
plot!(F,f,plot_ext=true,subplot=2)
Plots.savefig(P,"images/highprecision.png")
```

![](images/highprecision.png)

![](images/highprecision.png)

# Frame Approximations in 2D

In higher dimensions, a basis can be any tensorproduct of (scaled) lower dimensional bases:


```julia
using StaticArrays
C = disk(1.0)\disk(0.3,SVector(0.2, 0.5))
B = FourierBasis(31,-1.3,1.3) ⊗ FourierBasis(31,-1.3,1.3)
f = (x,y)->exp(x+y)
F = Fun(f,B,C)

P = heatmap(F,plot_ext=true,layout=2,aspect_ratio=1)
plot!(F,f,plot_ext=true,subplot=2,aspect_ratio=1)
Plots.savefig(P,"images/deathstar.png")
```

    Tuple{Int64, Int64}


![](images/deathstar.png)

Even fractal domains are not a problem:


```julia
B = FourierBasis(31,-1.0,0.35) ⊗ FourierBasis(31,-0.65,0.65)
f = (x,y)->cos(10*x*y)
F = Fun(f, B, mandelbrot())

P = heatmap(F,plot_ext=true,layout=2,aspect_ratio=1)
plot!(F,f,plot_ext=true,aspect_ratio=1,subplot=2)
Plots.savefig(P,"images/mandelbrot")
```

![](images/mandelbrot.png)
