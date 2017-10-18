
using BasisFunctions
using FrameFun
using Plots
using Domains
using Interact, Reactive

D = interval(-0.5,0.5)
B = FourierBasis(25)
F = Fun(x->x,B,D)
USV = svd(matrix(F))
plot(USV[2])

@manipulate for i=1:length(USV[2])
    F = SetFun(extensionframe(D,B),USV[3][:,i])
    plot(F,ylim=[-4,4],plot_ext=true)
end

D2 = mandelbrot()
B2 = FourierBasis(15,-1.0,0.35)⊗FourierBasis(15,-0.65,0.65)
F2 = Fun((x,y)->x+y,B2,D2)

USV2=svd(matrix(F2))
plot(USV2[2])

n = length(USV2[2])
@manipulate for i=1:n
    F2 = SetFun(extensionframe(D2,B2),reshape(USV2[3][:,i],round(Int,sqrt(n)),round(Int,sqrt(n))))
    plot(F2,n=101,title="l=$(USV2[2][i])",plot_ext=true)
end




