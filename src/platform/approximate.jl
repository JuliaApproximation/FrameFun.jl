
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
function approximate(::DictionaryOperatorStyle, samplingstyle::SamplingStyle, solverstyle::SolverStyle, fun, ap;
            verbose = false, options...)
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
    C = solve(solverstyle, ap, A, B; S=S, verbose=verbose, samplingstyle=samplingstyle, options...)
    F = DictFun(dictionary(ap), C)
    res = norm(A*C-B)
    verbose && println("Approximate: ended with residual $res\n")
    F, A, B, C, S, L
end

# Construct the approximation problem and solve it
function approximate(pstyle::GenericOperatorStyle, samplingstyle::SamplingStyle, solverstyle::SolverStyle, fun, ap;
            verbose = false, options...)
    dict = dictionary(ap)

    verbose && println("Approximate: using problemstyle style $pstyle")
    verbose && println("Approximate: using sampling style $samplingstyle")
	verbose && println("Approximate: using solver style $solverstyle")
    # Trigger computation of sampling parameter L first
    L = samplingparameter(samplingstyle, ap; verbose=verbose, options...)
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
    F = DictFun(dictionary(ap), C)
    res = norm(apply(S,A; options...)*C-B)
    verbose && println("Approximate: ended with residual $res\n")
    F, A, B, C, S, L
end


solve(style::SolverStyle, ap, A::DictionaryOperator, B; options...) =
    solver(style, ap, A; B=B, options...) * B
solve(style::SolverStyle, ap, A::AbstractOperator, fun; options...) =
    apply(solver(style, ap, A; options...), fun; options...)

solver(::InverseStyle, ap, A; options...) = inv(A)
solver(::DirectStyle, ap, A; options...) = directsolver(A; options...)

solver(style::AZStyle, ap, A; problemstyle=ProblemStyle(ap), options...) =
    solver(problemstyle, style, ap, A; options...)

solver(::DictionaryOperatorStyle, style::AZStyle, ap, A; options...) =
    solver(style, ap, A, AZ_Zt(ap; options...); options...)

function solver(::AZStyle, ap, A, Zt;
            B=nothing, smallcoefficients=false, smallcoefficients_atol=NaN, smallcoefficients_rtol=NaN, verbose=false, options...)
    if smallcoefficients
        Q = discrete_normalization(ap; options...)
        normF = abs(sqrt(sum(Q * B.^2)))
        if !isnan(smallcoefficients_rtol)
            verbose && println("Change smallcoefficients relative tolerance to absolute tolerance rtol*||f||")
            smallcoefficients_atol = smallcoefficients_rtol*normF
            smallcoefficients_rtol = NaN
        end

        AZSolver(A, Zt; smallcoefficients=smallcoefficients, smallcoefficients_rtol=smallcoefficients_rtol,
                smallcoefficients_atol=smallcoefficients_atol, verbose=verbose, options...)
    else
        AZSolver(A, Zt; B=B, verbose=verbose, options...)
    end
end

solver(style::AZSmoothStyle, ap, A; options...) = solver(style, ap, A, AZ_Zt(ap; options...); options...)
solver(::AZSmoothStyle, ap, A, Zt; options...) = AZSolver_with_smoothing(A, Zt; options...)

solver(::DualStyle, ap, A; options...) = dualdiscretization(ap; options...)

solver(solverstyle::ProductSolverStyle, ap, A; samplingstyle=SamplingStyle(ap), options...) =
    solver(solverstyle, samplingstyle, ap, A; options...)
solver(solverstyle::ProductSolverStyle, samplingstyle::ProductSamplingStyle, ap, A; S, options...) =
    TensorProductOperator(
		map( (ap_el,Ael,Sel,style,sstyle) -> solver(style, ap_el, Ael; S=Sel, samplingstyle=sstyle, options...),
				elements(ap), productelements(A), productelements(S), solverstyle.styles, samplingstyle.styles)...
	)

solver(::TridiagonalProlateStyle, ap, A; scaling = Zt_scaling_factor(dictionary(ap), A), options...) =
    FE_TridiagonalSolver(A, scaling; options...)

solver(pstyle::GenericOperatorStyle, solverstyle::AZStyle, ap::ApproximationProblem, A; options...) =
    GenericAZSolver(ap, A; options...)
