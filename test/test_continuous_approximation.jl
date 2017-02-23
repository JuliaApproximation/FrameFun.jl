using BasisFunctions
using FrameFun


for basistype in (ChebyshevBasis, FourierBasis, PeriodicBSplineBasis), T in (Float32, Float64)
  tol = sqrt(eps(T))

  B = instantiate(basistype,11, T)
  println(B)
  D = Interval(left(B),right(B))
  frame = extensionframe(B,D)
  f = x->B[1](x)

  frameop = approximation_operator(frame; discrete=false, abstol=tol)
  framopmatrix = matrix(frameop)
  basisop = approximation_operator(B; discrete=false, abstol=tol)
  basisopmatrix = matrix(basisop)
  norm(matrix(frameop-basisop))
  println(norm(matrix(frameop-basisop)))

  framesol = *(frameop,f; discrete=false, abstol=tol)

  basissol = *(basisop,f; discrete=false, abstol=tol)

  norm(framesol-basissol)
  println(norm(framesol-basissol))

  frameF = approximate(frame, f; discrete=false, abstol=tol)
  frameFcoef = coefficients(frameF)
  basisF = approximate(B, f; discrete=false, abstol=tol)
  basisFcoef = coefficients(basisF)

  norm(frameFcoef-basisFcoef)
  println(norm(frameFcoef-basisFcoef))
end

using Plots

plot(frameF)
plot!(basisF, linestyle=:dot, ylims=[-1,10])
