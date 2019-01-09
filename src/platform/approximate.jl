
###################
# Helper functions
###################

determine_return_type(fun, ::Type{S}) where {S} = Base.Core.Compiler.return_type(fun, (S,))
determine_return_type(fun, S::Type{<:Tuple}) = Base.Core.Compiler.return_type(fun, S)

function promote_dictionary(dict, fun)
    T = promote_type(codomaintype(dict), determine_return_type(fun, domaintype(dict)))
    promote_coefficienttype(dict, T)
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

SamplingStyle(ap::ApproximationProblem, M::Nothing) = SamplingStyle(ap)
SamplingStyle(ap::ApproximationProblem, M) = OversamplingStyle()

function approximate(fun, ap::ApproximationProblem;
        M = nothing,
        samplingstyle = SamplingStyle(ap, M),
        solverstyle = SolverStyle(ap, samplingstyle), options...)
    if M == nothing
        approximate(samplingstyle, solverstyle, fun, ap; options...)
    else
        approximate(samplingstyle, solverstyle, fun, ap; M=M, options...)
    end
end


# Construct the approximation problem and solve it
function approximate(samplingstyle::SamplingStyle, solverstyle::SolverStyle, fun, ap; verbose = false, options...)
    dict = dictionary(ap)

    S = samplingoperator(samplingstyle, ap; verbose=verbose, options...)
    if verbose
        println("Approximate: using sampling style $samplingstyle")
        showsamplinginformation(samplingstyle, dict, S)
    end

    A = apply(S, dict)
    B = apply(S, fun)
    verbose && println("Approximate: using solver style $solverstyle")
    C = solve(solverstyle, ap, A, B; S=S, verbose=verbose, options...)
    F = DictFun(dictionary(ap), C)
    err = norm(A*C-B)
    verbose && println("Approximate: ended with residual $err\n")
    F, A, B, C, S
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
    println(BasisFunctions.print_strings(("Fun: discrete approximation with grid:", BasisFunctions.strings(g))))
end



solve(style::SolverStyle, ap, A, B; options...) =
    solver(style, ap, A; options...) * B

solver(::InverseStyle, ap, A; options...) = inv(A)
solver(::DirectStyle, ap, A; options...) = directsolver(A; options...)

function solver(::AZStyle, ap, A;
            S = nothing,
            Zt = AZ_Zt(ap, S),
            options...)
    AZSolver(A, Zt; options...)
end

function solver(::AZSmoothStyle, ap, A;
            S = nothing,
            Zt = AZ_Zt(ap, S),
            options...)
    # AZSmoothSolver(A, Zt; options...)
    AZSolver_with_smoothing(A, Zt; options...)
end


function solver(::DualStyle, ap, A; S, options...)
    Stilde = dualsamplingoperator(ap, S)
    dtilde = dualdictionary(ap)
    B = apply(Stilde, dtilde)
    B'
end

solver(::TridiagonalProlateStyle, ap, A; scaling = Zt_scaling_factor(dictionary(ap), A), options...) =
    FE_TridiagonalSolver(A, scaling; options...)
