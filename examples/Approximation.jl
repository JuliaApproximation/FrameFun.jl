
using BasisFunctions
using FrameFun
using Domains
using Plots;gr()
using StaticArrays

B = FourierBasis(61,-1,1)
D = interval(-0.5,0.5)
f1 = x->cos(3*x)
F1 = Fun(f1,B,D)

F1(0.5)

f1(0.5)

plot(F1,layout=2)
plot!(F1,f1,subplot=2)



B2 = ChebyshevBasis(130)
D = interval()/2
f2 = x->cos(80*x)
F2 = Fun(f2, B2, D)

F2(0.1)

f2(0.1)

plot(F2, layout=2, plot_ext=true)
plot!(F2,f2, subplot=2)



f3 = x->cos(10*x.^2)
B = FourierBasis(41,-1,1)
D = interval(-1.0,-0.5)∪interval(-0.2,0.5)
F3 = Fun(f3,B,D)

l = @layout [Plots.grid(1,1); Plots.grid(1,2)]
plot(F3, layout=l)
plot!(F3, subplot=2, plot_ext=true)
plot!(F3,f3, subplot=3)



B = FourierBasis(61,BigFloat)
D = interval(0.,0.5)
fh = x->x
Fh = Fun(fh,B,D)

pt = 3//10
abs(Fh(pt)-fh(pt))

plot(Fh, layout=2)
plot!(Fh,fh, subplot=2)



D = interval(0.,.5)^2

B = FourierBasis(100)⊗FourierBasis(100)
f = (x,y)->exp(y*2*x)
F = Fun(f,B,D)

plot(F)

plot(F,f)

C = disk() \ disk(0.3,SVector(0.2, 0.5))

plot(C)

B = FourierBasis(41,-1.3,1.3) ⊗ FourierBasis(41,-1.3,1.3)
fC(x,y) = exp(y+x)
F = Fun(fC,B,C)

F(0,0.4)

fC(0, 0.4)

plot(F)

heatmap(F,fC)



D=FrameFun.mandelbrot()

plot(D)

B = FourierBasis(31,-1.0,0.35) ⊗ FourierBasis(31,-0.65,0.65)
fm(x,y) = cos(10*x*y)
F = Fun(fm, B, D)

heatmap(F)

contourf(F,fm)




