
###################
# Helper functions
###################

determine_return_type(fun, S) = Base.promote_op(fun, S)
determine_return_type(fun, S::Type{Tuple{A}}) where {A} = Base.promote_op(fun, A)
determine_return_type(fun, S::Type{Tuple{A,B}}) where {A,B} = Base.promote_op(fun, A, B)
determine_return_type(fun, S::Type{Tuple{A,B,C}}) where {A,B,C} = Base.promote_op(fun, A, B, C)
determine_return_type(fun, S::Type{SVector{1,A}}) where {A} = Base.promote_op(fun, A)
determine_return_type(fun, S::Type{SVector{2,A}}) where {A} = Base.promote_op(fun, A, A)
determine_return_type(fun, S::Type{SVector{3,A}}) where {A} = Base.promote_op(fun, A, A, A)
determine_return_type(fun, S::Type{SVector{4,A}}) where {A} = Base.promote_op(fun, A, A, A, A)


promote_dictionary(dict, fun) = _promote_dictionary(dict, fun, determine_return_type(fun, domaintype(dict)))
_promote_dictionary(dict, fun, ::Union{}) = dict
_promote_dictionary(dict, fun, ::Type{T}) where {T} = BasisFunctions.ensure_coefficienttype(promote_type(codomaintype(dict),T), dict)


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
	elseif directsolver == :regsvd
        verbose && println("Using direct solver: regularized SVD")
        regularized_SVD_solver(A; verbose=verbose, options...)
    else
        DS = directsolver(A; verbose=verbose, options...)
        verbose && println("Using direct solver: $DS")
        DS
    end
end

function iterativesolver(A; iterativesolver = :lsqr, verbose = false, options...)
    if iterativesolver == :lsqr
        verbose && println("Using iterative solver: LSQR")
        LSQR_solver(A; verbose=verbose, options...)
    elseif iterativesolver == :lsmr
        verbose && println("Using iterative solver: LSMR")
        LSMR_solver(A; verbose=verbose, options...)
    else
        DS = iterativesolver(A)
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

guess_coefficienttype(dict, fun) = promote_type(codomaintype(dict), determine_return_type(fun, domaintype(dict)))

function Fun(fun, dict::Dictionary, args...;
        coefficienttype = guess_coefficienttype(dict, fun),
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
        show(dictionary(ap))
		println()
    end
    F, A, B, C, S = approximate(fun, ap; verbose=verbose, options...)
    F
end

function Fun(fun, ap::AdaptiveApproximation; verbose = false, options...)
    verbose && println("Fun: adaptive approximation using platform $(ap.platform)")
    approximate(fun, ap; verbose=verbose, options...)
end


function Fun(dict::Dictionary, coefficients::AbstractArray; verbose = false, options...)
	verbose && println("Fun: coefficients are given, no approximation problem is solved")
	Expansion(dict, coefficients)
end


############################
# The approximate function
############################


approximate(fun, dict::Dictionary, args...; options...) =
    approximate(fun, approximationproblem(dict, args...); options...)
approximate(fun, platform::Platform, args...; options...) =
    approximate(fun, approximationproblem(platform, args...); options...)

approximate(fun, ap::ApproximationProblem;
            samplingstyle = SamplingStyle(ap),
            solverstyle = SolverStyle(samplingstyle, ap),
            problemstyle = ProblemStyle(ap), options...) =
    approximate(problemstyle, promote(samplingstyle, solverstyle)..., fun, ap; options...)

# If the sampling has product structure and a single solverstyle is specified,
# we turn the solverstyle into a product style as well.
promote(samplingstyle::SamplingStyle, solverstyle::SolverStyle) = (samplingstyle,solverstyle)
promote(samplingstyle::ProductSamplingStyle, solverstyle::ProductSolverStyle) = (samplingstyle,solverstyle)
promote(samplingstyle::ProductSamplingStyle, solverstyle::SolverStyle) =
	(samplingstyle,ProductSolverStyle(map(x->solverstyle, samplingstyle.styles)))

# Construct the approximation problem and solve it
function approximate(::DictionaryOperatorStyle, samplingstyle::SamplingStyle, solverstyle::SolverStyle, fun, ap::ApproximationProblem;
            verbose = false, options...)
    dict = dictionary(ap)

    verbose && println("Approximate: using sampling style $samplingstyle")
	verbose && println("Approximate: using solver style $solverstyle")

    S = samplingoperator(samplingstyle, ap; verbose=verbose, options...)
	if verbose
		println()
		println("Approximate: sampling operator with size $(size(S)) is")
		println(S)
		println()
	end
    A, B = normalized_discretization(fun, samplingstyle, ap, S; verbose=verbose, options...)
    C = solve(solverstyle, ap, A, B; S=S, verbose=verbose, samplingstyle=samplingstyle, options...)
    F = Expansion(dictionary(ap), C)
    res = norm(A*C-B)
    verbose && println("Approximate: ended with residual $res\n")
    F, A, B, C, S, samplingparameter(samplingstyle, ap; verbose=verbose, options...)
end

# Construct the approximation problem and solve it
function approximate(pstyle::GenericOperatorStyle, samplingstyle::SamplingStyle, solverstyle::SolverStyle, fun, ap;
            verbose = false, options...)
    dict = dictionary(ap)

    verbose && println("Approximate: using problemstyle style $pstyle")
    verbose && println("Approximate: using sampling style $samplingstyle")
	verbose && println("Approximate: using solver style $solverstyle")

    S = samplingoperator(pstyle, ap; samplingstyle=samplingstyle, verbose=verbose, options...)
	if verbose
		println()
		println("Approximate: sampling operator with size $(size(S)) is")
		println(S)
		println()
	end

    A = AZ_A(pstyle, ap; samplingstyle=samplingstyle, options...)
    B = apply(S, fun; options...) # calculated just to keep the same interface.
    C = solve(solverstyle, ap, A, fun; verbose=verbose, problemstyle=pstyle, samplingstyle=samplingstyle, options...)
    F = Expansion(dictionary(ap), C)
    res = norm(apply(S,A; options...)*C-B)
    verbose && println("Approximate: ended with residual $res\n")
    F, A, B, C, S, samplingparameter(samplingstyle, ap; verbose=verbose, options...)
end


solve(style::SolverStyle, ap, A::DictionaryOperator, B; options...) =
    solver(style, ap, A; B=B, options...) * B
solve(style::SolverStyle, ap, A::AbstractOperator, fun; options...) =
    apply(solver(style, ap, A; options...), fun; options...)
