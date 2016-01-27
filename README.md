[![Build Status](https://travis-ci.org/daanhb/FrameFuns.jl.svg?branch=master)](https://travis-ci.org/daanhb/FrameFuns.jl)

FrameFuns
=========

Exploring practical possibilities of approximating functions with frames rather than with a basis. The package is heavily inspired by the Chebfun project and the Julia package ApproxFun.

# Frame Approximations in 1D

After choosing a suitable Basis and Domain, any function can be approximated in the resulting frame:
```julia
using FrameFuns
B = FourierBasis(61)
D = Interval(-0.5,0.5)
f(x) = x
F = ExpFun(f,D)

plot_expansion(F); plot_error(F,f)
```

![](images/lowprecision.png)

The bases support any AbstractFloat subtype, so high precision approximations are straightforward:

```julia
B = FourierBasis(61,BigFloat)
F = ExpFun(f,D)

plot_expansion(F); plot_error(F,f)
```

![](images/highprecision.png)

# Frame Approximations in 2D

In higher dimensions, a basis can be any tensorproduct of (scaled) lower dimensional bases:
```julia
C = Disk(1.0)- Disk(0.3,[0.2; 0.5])
B = rescale(FourierBasis(31),-1.3,1.3) ⊗ rescale(FourierBasis(31),-1.3,1.3)
f(x,y) = exp(x+y)
F = Fun(f,B,C)

plot_image(F); plot_error(F,f)
```

![](images/deathstar.png)

Basis types can easily be mixed and matched, and domains can be arbitrarily complex:

```julia
dom = randomcircles(10)
B = FourierBasis(31) ⊗ ChebyshevBasis(31)
f(x,y) = cos(20*x+22*y)
F = Fun(f,B,dom)

plot_image(F), plot_error(F,f)
```

![](images/circles.png)

Even fractal domains are not a problem:

```julia
B = rescale(FourierBasis(31),-1.0,0.35) ⊗ rescale(FourierBasis(31),-0.65,0.65)
f(x,y) = cos(10*x*y)
F = Fun(f, B, Mandelbrot())

plot_image(F), plot_error(F,f)
```

![](images/mandelbrot.png)
