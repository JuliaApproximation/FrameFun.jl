
using BasisFunctions
using FrameFun
using Plots
using DomainSets

B = ChebyshevBasis(10,-2.0,2.0)
D = interval(-1.0,1.0)
f = x->x^2
F = Fun(f,B,D)

plot(∫(F),layout=4)
plot!(F',subplot=2)
plot!(F,subplot=3)
plot!(∫(F'),subplot=4)

B = FourierBasis(100,-1,1)⊗FourierBasis(100,-1,1)
D = interval(-0.5,0.5)×interval(-0.5,0.5)
f = (x,y)->x*y
F = Fun(f,B,D)

heatmap(F,layout=2)
heatmap!(∂x(F),subplot=2)

plot(F,f,layout=2)
plot!(∂x(F),(x,y)->y,subplot=2)
