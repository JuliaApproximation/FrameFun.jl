
using BasisFunctions
using FrameFun
using Plots
using DomainSets

FE = FrameFun
BA = BasisFunctions

B = Fourier(41,-1,1)
Dom = interval(-0.5,0.5)

diff = IdentityOperator(B)
BC = DirichletBC(x->0);

fD = x->x;
Diff = differentiation_operator(B)*differentiation_operator(B)
DE = DiffEquation(B,Dom,Diff,fD, (BC,BC));

FD = solve(DE)

# Exact solution
solD = x->x^3/6 - x/24;

plot(FD,layout=4,title="Solution")
plot!(FD'',subplot=2,title="Second derivative")
plot!(FD,solD,subplot=3,title="Solution error")
plot!(FD'',fD,subplot=4,title="Derivative error")

diff = differentiation_operator(B)
NeumannBC(x->0.)

fN = x->x;
Diff = differentiation_operator(B)*differentiation_operator(B)
DE = DiffEquation(B,Dom,Diff,fN, (BC,));

FN = solve(DE)

# Exact solution
solN = x->x^3/6-x/8

plot(FN,layout=4,title="Solution")
plot!(FN'',subplot=2,title="Second derivative")
plot!(FN,solN,subplot=3,title="Solution error")
plot!(FN'',fN,subplot=4,title="Derivative error")

B2 = Fourier(31,-1.0,1.0)⊗Fourier(31,-1.0,1.0)
D2 = disk(0.8)\disk(0.3)\cube([-0.15,-1.0],[0.15,0.0])

diff2 = IdentityOperator(B2)
df2D = (x,y) -> x-y;
BC2 = DirichletBC(df2D,euclideanspace(Val{2}()))

f2D = (x,y)->0;
Diff2 = differentiation_operator(B2,(2,0))+differentiation_operator(B2,(0,2))
DE2 = DiffEquation(B2,D2,Diff2,f2D, (BC2,));

F2D = solve(DE2,solver=DirectSolver)  

heatmap(F2D)

B2 = Fourier(21,-1.0,1.0)⊗Fourier(21,-1.0,1.0)
D2 = cube([-1.0,-0.5],[1.0,0.5])

diff2 = differentiation_operator(B2,(0,1))
df2N = (x,y)->0;
BC2 = NeumannBC(df2N,euclideanspace(Val{2}()))

f2N = (x,y)->cos(2*pi*(x+y));
Diff2 = differentiation_operator(B2,(2,0))+differentiation_operator(B2,(0,2))
DE2 = DiffEquation(B2,D2,Diff2,f2N, (BC2,));

F2N = solve(DE2,solver=AZSolver)

heatmap(F2N)

using StaticArrays
B2 = Fourier(41,-1.0,1.0)⊗Fourier(41,-1.0,1.0)
D2 = disk(0.75)\disk(0.2,SVector(0.3,-0.3))

diff2 = IdentityOperator(B2)
df2H = (x,y)->0;
BC = NeumannBC(df2H,euclideanspace(Val{2}()))

f2H = (x,y)->exp(-200*((x+0.3)^2+(y-0.3)^2))
Diff2 = differentiation_operator(B2,(2,0))+differentiation_operator(B2,(0,2))+1000*IdentityOperator(B2)
DE2 = DiffEquation(B2,D2,Diff2,f2H, (BC2,));

F2H = FrameFun.solve(DE2,solver=AZSolver)

heatmap(F2H)
