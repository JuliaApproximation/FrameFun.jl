
import BasisFunctions: approximate
export Fun

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


guess_coefficienttype(dict, fun) =
	_guess_coefficienttype(codomaintype(dict), determine_return_type(fun, domaintype(dict)))

_guess_coefficienttype(::Type{Any}, ::Type{Any}) = Any
_guess_coefficienttype(::Type{T}, ::Type{Any}) where {T} = T
_guess_coefficienttype(::Type{Any}, ::Type{T}) where {T} = Any
_guess_coefficienttype(::Type{S}, ::Type{T}) where {S,T} = promote_type(S,T)

"""
`Fun` is used to approximate functions in the same way as `approximate` with
the same interface. It discards all outputs of `approximate` and only
returns the function approximation.
"""
function Fun(fun, dict::Dictionary, args...;
        coefficienttype = guess_coefficienttype(dict, fun),
        options...)
	dict2 = BasisFunctions.ensure_coefficienttype(coefficienttype, dict)
    ap = approximationproblem(fun, dict2, args...; options...)
    Fun(ap; options...)
end

function Fun(fun, platform::Platform, args...; options...)
    ap = approximationproblem(fun, platform, args...)
    Fun(ap; options...)
end

Fun(fun, ap::ApproximationProblem; options...) = Fun(_ap(fun, ap); options...)

function Fun(ap::ApproximationProblem; verbose = false, options...)
    if verbose
        println("Fun: using the following dictionary:")
        show(dictionary(ap))
		println()
    end
    F, A, B, C, S = approximate(ap; verbose=verbose, options...)
    F
end

function Fun(ap::AdaptiveApproximation; verbose = false, options...)
    verbose && println("Fun: adaptive approximation using platform $(platform(ap))")
    approximate(ap; verbose=verbose, options...)
end


function Fun(dict::Dictionary, coefficients::AbstractArray; verbose = false, options...)
	verbose && println("Fun: coefficients are given, no approximation problem is solved")
	Expansion(dict, coefficients)
end


############################
# The approximate function
############################


approximate(fun, dict::Dictionary, args...; options...) =
    approximate(approximationproblem(fun, dict, args...; options...); options...)
approximate(fun, platform::Platform, args...; options...) =
    approximate(approximationproblem(fun, platform, args...; options...); options...)

approximate(ap::ApproximationProblem;
            samplingstyle = SamplingStyle(ap),
            solverstyle = SolverStyle(ap), options...) =
    approximate(promote(samplingstyle, solverstyle)..., ap; options...)

# If the sampling has product structure and a single solverstyle is specified,
# we turn the solverstyle into a product style as well.
promote(samplingstyle::SamplingStyle, solverstyle::SolverStyle) = (samplingstyle,solverstyle)
promote(samplingstyle::ProductSamplingStyle, solverstyle::ProductSolverStyle) = (samplingstyle,solverstyle)
promote(samplingstyle::ProductSamplingStyle, solverstyle::SolverStyle) =
	(samplingstyle,ProductSolverStyle(map(x->solverstyle, samplingstyle.styles)))

# Construct the approximation problem and solve it
function approximate(samplingstyle::SamplingStyle, solverstyle::SolverStyle, ap::ApproximationProblem;
            verbose = false, options...)
    dict = dictionary(ap)

    verbose && println("Approximate: using sampling style $samplingstyle")
	verbose && println("Approximate: using solver style $solverstyle")

    S = sampling_operator(samplingstyle, ap; verbose=verbose, options...)
	if verbose
		println()
		println("Approximate: sampling operator with size $(size(S)) is")
		println(S)
		println()
	end
	A = discretization(ap; options...)
	b = sample_data(ap; options...)
    coef = solve(solverstyle, ap, A, b; S=S, verbose=verbose, samplingstyle=samplingstyle, options...)
    F = Expansion(dictionary(ap), coef)
    res = norm(A*coef-b)
    verbose && println("Approximate: ended with residual $res\n")
    F, A, b, coef, S, samplingparameter(ap)
end

export solve
solve(style::SolverStyle, ap, A::DictionaryOperator, B; options...) =
    solver(style, ap, A; B=B, options...) * B
solve(style::SolverStyle, ap, A::AbstractOperator, fun; options...) =
    apply(solver(style, ap, A; options...), fun; options...)
