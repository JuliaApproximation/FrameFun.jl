module WeightedApproximation

using BasisFunctions, LinearAlgebra
using FrameFun

import BasisFunctions: approximate

struct WeightedApproximationStyle <: ProblemStyle
    weight
end

function approximate(pstyle::WeightedApproximationStyle, samplingstyle::SamplingStyle, solverstyle::SolverStyle, fun, ap;
            verbose = false, options...)
    dict = dictionary(ap)

    verbose && println("Approximate: using problem style $pstyle")
    verbose && println("Approximate: using sampling style $samplingstyle")
	verbose && println("Approximate: using solver style $solverstyle")
    # Trigger computation of sampling parameter L first
    L = samplingparameter(samplingstyle, ap; verbose=verbose, options...)
    S = samplingoperator(samplingstyle, ap; verbose=verbose, options...)
	if verbose
		println()
		println("Approximate: sampling operator with size $(size(S)) is")
		println(S)
		println()
	end
    A = apply(S, dict; options...)
    B = apply(S, fun; options...)
    W = DiagonalOperator(apply(S, pstyle.weight; options...))

    C = solve(solverstyle, ap, W*A, W*B; problemstyle=pstyle, S=S, verbose=verbose, samplingstyle=samplingstyle, options...)
    F = Expansion(dictionary(ap), C)
    res = norm(W*(A*C-B))
    verbose && println("Approximate: ended with residual $res\n")
    F, A, B, C, S, L
end


export WFun
"""
Creates an approximate solution to `fun` using a weighed L2 measure, given by the `weight`.
"""
WFun(fun, weight, args...; options...) =
    Fun(fun, args...;
        problemstyle=WeightedApproximationStyle(weight), solverstyle=AZStyle(), samplingstyle=OversamplingStyle(),
        options...)

solver(pstyle::WeightedApproximationStyle, style::AZStyle, ap, A; options...) =
    solver(style, ap, A, AZ_Zt(pstyle, ap; options...); options...)

AZ_A(pstyle::WeightedApproximationStyle, ap; options...) =
    WAZ_W(pstyle, ap; options...)*AZ_A(DictionaryOperatorStyle(), ap; options...)

AZ_Z(pstyle::WeightedApproximationStyle, ap; weight_epsilon=1e-12, options...) =
    pinv(WAZ_W(pstyle, ap; options...), weight_epsilon)*AZ_Z(DictionaryOperatorStyle(), ap; options...)

WAZ_W(pstyle, ap::ApproximationProblem; samplingstyle=SamplingStyle(ap), options...) =
    DiagonalOperator(apply(samplingoperator(samplingstyle, ap; options...), pstyle.weight; options...))

end
