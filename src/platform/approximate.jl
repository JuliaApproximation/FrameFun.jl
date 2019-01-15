
###################
# Helper functions
###################

determine_return_type(fun, ::Type{S}) where {S} = Base.Core.Compiler.return_type(fun, (S,))
determine_return_type(fun, S::Type{<:Tuple}) = Base.Core.Compiler.return_type(fun, S)

function promote_dictionary(dict, fun)
    T = promote_type(codomaintype(dict), determine_return_type(fun, domaintype(dict)))
    promote_coefficienttype(dict, T)
end

function directsolver(A::DiagonalOperator; verbose = false, options...)
    verbose && println("Using direct solver: inverse of diagonal matrix")
    inv(A)
end

function directsolver(A; directsolver = :svd, verbose = false, options...)
    if directsolver == :svd
        verbose && println("Using direct solver: SVD")
        SVD_solver(A)
    elseif directsolver == :qr
        verbose && println("Using direct solver: QR")
        QR_solver(A)
    else
        DS = directsolver(A)
        verbose && println("Using direct solver: $DS")
        DS
    end
end



####################
# The Fun interface
####################

# The `Fun` interface first turns the approximation problem into a subtype of
# `ApproximationProblem`. Next, it calls `approximate`. Finally, it discards all
# the operators that were computed and simply returns the function.
# Users wanting to access the operators can call `approximate` directly.

function Fun(fun, dict::Dictionary, args...;
        coefficienttype = promote_type(codomaintype(dict), determine_return_type(fun, domaintype(dict))),
        options...)
    ap = approximationproblem(coefficienttype, dict, args...)
    Fun(fun, ap; options...)
end

function Fun(fun, platform::Platform, args...; options...)
    ap = approximationproblem(platform, args...)
    Fun(fun, ap; options...)
end

function Fun(fun, ap::ApproximationProblem; verbose = false, options...)
    if verbose
        println("Fun: using the following dictionary:")
        println("$(dictionary(ap))")
    end
    F, A, B, C, S = approximate(fun, ap; verbose=verbose, options...)
    F
end

function Fun(fun, ap::AdaptiveApproximation; verbose = false, options...)
    verbose && println("Fun: adaptive approximation using platform $(ap.platform)")
    approximate(fun, ap; verbose=verbose, options...)
end


Fun(dict::Dictionary, coefficients::AbstractArray) = DictFun(dict, coefficients)


############################
# The approximate function
############################


approximate(fun, dict::Dictionary, args...; options...) =
    approximate(fun, approximationproblem(dict, args...); options...)
approximate(fun, platform::Platform, args...; options...) =
    approximate(fun, approximationproblem(platform, args...); options...)

approximate(fun, ap::ApproximationProblem;
            samplingstyle = SamplingStyle(ap),
            solverstyle = SolverStyle(ap, samplingstyle), options...) =
    approximate(samplingstyle, solverstyle, fun, ap; options...)


# Construct the approximation problem and solve it
function approximate(samplingstyle::SamplingStyle, solverstyle::SolverStyle, fun, ap; verbose = false, options...)
    dict = dictionary(ap)

    verbose && println("Approximate: using sampling style $samplingstyle")
    # Trigger computation of sampling parameter L first
    L = samplingparameter(samplingstyle, ap; verbose=verbose, options...)
    S = samplingoperator(samplingstyle, ap; verbose=verbose, options...)
    verbose && showsamplinginformation(samplingstyle, dict, S)

    A = apply(S, dict)
    B = apply(S, fun)
    verbose && println("Approximate: using solver style $solverstyle")
    C = solve(solverstyle, ap, A, B; S=S, verbose=verbose, options...)
    F = DictFun(dictionary(ap), C)
    err = norm(A*C-B)
    verbose && println("Approximate: ended with residual $err\n")
    F, A, B, C, S, L
end


showsamplinginformation(dstyle::DiscreteStyle, dict::Dictionary, S::AbstractOperator) =
    showsamplinginformation(dstyle, element(S, numelements(S)))

function showsamplinginformation(::DiscreteStyle, dict::Dictionary, S::GridSampling)
    g = grid(dest(S))
    N = length(dict)
    M = length(g)
    if M > N
        println("Fun: oversampling with N=$N and M=$M")
    end
    println(BasisFunctions.print_strings(("Fun: discrete sampling operator:", BasisFunctions.strings(S))))
end


function showsamplinginformation(dstyle::GramStyle, dict::Dictionary, S::ProjectionSampling)
    println(BasisFunctions.print_strings(("Fun: continuous approximation with projection:", BasisFunctions.strings(S))))
end


solve(style::SolverStyle, ap, A, B; options...) =
    solver(style, ap, A; options...) * B

solver(::InverseStyle, ap, A; options...) = inv(A)
solver(::DirectStyle, ap, A; options...) = directsolver(A; options...)

solver(style::AZStyle, ap, A; options...) = solver(style, ap, A, AZ_Zt(ap; options...); options...)
solver(::AZStyle, ap, A, Zt; options...) = AZSolver(A, Zt; options...)

solver(style::AZSmoothStyle, ap, A; options...) = solver(style, ap, A, AZ_Zt(ap; options...); options...)
solver(::AZSmoothStyle, ap, A, Zt; options...) = AZSolver_with_smoothing(A, Zt; options...)

solver(::DualStyle, ap, A; options...) = dualdiscretization(ap; options...)

solver(::TridiagonalProlateStyle, ap, A; scaling = Zt_scaling_factor(dictionary(ap), A), options...) =
    FE_TridiagonalSolver(A, scaling; options...)
