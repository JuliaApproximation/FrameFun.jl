using BasisFunctions
using FrameFun
using Base.Test


for basistype in (ChebyshevBasis, FourierBasis, BSplineTranslatesBasis), T in (Float32, Float64,)
  tol = sqrt(eps(T))

  B = instantiate(basistype,11, T)
  # println(B)
  D = interval(left(B),right(B))
  frame = extensionframe(B,D)
  f = x->B[1](x)

  frameop = approximation_operator(frame; discrete=false, abstol=tol)
  # framopmatrix = matrix(frameop)
  basisop = approximation_operator(B; discrete=false, abstol=tol)
  # basisopmatrix = matrix(basisop)

  framesol = *(frameop,f; discrete=false, abstol=tol)

  basissol = *(basisop,f; discrete=false, abstol=tol)

  @test norm(framesol-basissol) < 10*tol
  # println(norm(framesol-basissol))

  frameF = approximate(frame, f; discrete=false, abstol=tol)
  frameFcoef = coefficients(frameF)
  basisF = approximate(B, f; discrete=false, abstol=tol)
  basisFcoef = coefficients(basisF)

  @test norm(frameFcoef-basisFcoef) < 10*tol
  # println(norm(frameFcoef-basisFcoef))
end

using Plots

plot(frameF)
plot!(basisF, linestyle=:dot, ylims=[-1,10])
