
using BasisFunctions
BA = BasisFunctions
using FrameFun
FE = FrameFun
using DomainSets
using Plots
using LinearAlgebra

gr()

function myFun(f::Function,P::Platform,i)
    AZS = AZSolver(A(P,i),Zt(P,i))
    F = DictFun(primal(P,i),AZS*f)
end

P = BA.fourier_platform()


f = x->cos(4*pi*x)
F = myFun(f,P,8)
plot(F)

f = x->exp(x)
F = myFun(f,P,8)
plot(F)

plot(svdvals(matrix(A(P,5)*Zt(P,5))),ylims=(-0.1,1.1))

D = interval(0,0.5)
P2a = BasisFunctions.fourier_platform(oversampling=4)
P2 = FrameFun.extension_frame_platform(P2a,D)

f = x->exp(x)
F2 = myFun(f,P2,8)
plot(F2)

plot(svdvals(matrix(A(P2,5)*Zt(P2,5))),ylims=(-0.1,1.1))

WSP = FrameFun.WeightedSumPlatform([x->sqrt(x),x->1],P2)

f = x->sqrt(x)*(1-x)-exp(x)
FSP = myFun(f,WSP,8)
plot(FSP,f)

FSPb = myFun(f,P2,8)
plot(FSP,f)
plot!(FSPb,f)

plot(svdvals(matrix(A(WSP,7)*Zt(WSP,7))),ylims=(-0.1,4.1))
