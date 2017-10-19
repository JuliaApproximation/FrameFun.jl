module test_suite
using BasisFunctions
using FrameFun
using Base.Test

function delimit(s::AbstractString)
    println()
    println("############")
    println("# ",s)
    println("############")
end

delimit("Continuous solver")

frameF = 0
basisF = 0
function continuous_solver_test()
    @testset for basistype in (ChebyshevBasis, FourierBasis, BSplineTranslatesBasis), T in (Float32, Float64,), solver in (FrameFun.FE_DirectSolver, FrameFun.ExactTruncatedSvdSolver)
        tol = sqrt(eps(T))

        B = instantiate(basistype,11, T)
        # println(B)
        D = interval(left(B),right(B))
        frame = extensionframe(B,D)
        f = x->B[1](x)

        frameop = approximation_operator(span(frame); discrete=false, solver=solver, abstol=tol)
        # framopmatrix = matrix(frameop)
        basisop = approximation_operator(span(B); discrete=false, solver=solver, abstol=tol)
        # basisopmatrix = matrix(basisop)

        framesol = *(frameop,f; discrete=false, abstol=tol)

        basissol = *(basisop,f; discrete=false, abstol=tol)

        @test norm(framesol-basissol) < 10*tol
        # println(norm(framesol-basissol))

        frameF = approximate(span(frame), f; discrete=false, solver=solver, abstol=tol)
        frameFcoef = coefficients(frameF)
        basisF = approximate(span(B), f; discrete=false, solver=solver, abstol=tol)
        basisFcoef = coefficients(basisF)

        @test norm(frameFcoef-basisFcoef) < 10*tol
        # println(norm(frameFcoef-basisFcoef))
    end
end
continuous_solver_test()

# using Plots
#
# plot(frameF)
# plot!(basisF, linestyle=:dot, linewidth=10, color=:red)
end
