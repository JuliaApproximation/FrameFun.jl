
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
	elseif directsolver == :regsvd
        verbose && println("Using direct solver: regularized SVD")
        regularized_SVD_solver(A; verbose=verbose, options...)
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
    approximate(promote(samplingstyle, solverstyle)..., fun, ap; options...)

# If the sampling has product structure and a single solverstyle is specified,
# we turn the solverstyle into a product style as well.
promote(samplingstyle::SamplingStyle, solverstyle::SolverStyle) = (samplingstyle,solverstyle)
promote(samplingstyle::ProductSamplingStyle, solverstyle::ProductSolverStyle) = (samplingstyle,solverstyle)
promote(samplingstyle::ProductSamplingStyle, solverstyle::SolverStyle) =
	(samplingstyle,ProductSolverStyle(map(x->solverstyle, samplingstyle.styles)))

# Construct the approximation problem and solve it
function approximate(samplingstyle::SamplingStyle, solverstyle::SolverStyle, fun, ap; verbose = false, options...)
    dict = dictionary(ap)

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
    C = solve(solverstyle, ap, A, B; S=S, verbose=verbose, options...)
    F = DictFun(dictionary(ap), C)
    res = norm(A*C-B)
    verbose && println("Approximate: ended with residual $res\n")
    F, A, B, C, S, L
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


solver(solverstyle::ProductSolverStyle, ap, A; S, options...) =
    TensorProductOperator(
		map( (ap_el,Ael,Sel,style) -> solver(style, ap_el, Ael; S=Sel, options...),
				elements(ap), productelements(A), productelements(S), solverstyle.styles)...
	)

solver(::TridiagonalProlateStyle, ap, A; scaling = Zt_scaling_factor(dictionary(ap), A), options...) =
    FE_TridiagonalSolver(A, scaling; options...)
