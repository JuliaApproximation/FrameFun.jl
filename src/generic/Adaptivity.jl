module Adaptivity
using Printf, ..Platforms, BasisFunctions, ..ApproximationProblems, ..DictFuns

import BasisFunctions: approximate

export ErrorStyle
abstract type ErrorStyle end
export ResidualStyle, OversampledResidual, RandomPoints, FNAStyle
struct ResidualStyle <: ErrorStyle end
struct OversampledResidual <: ErrorStyle end
struct RandomPoints <: ErrorStyle end
struct FNAStyle{CREL,EREL} <: ErrorStyle
    FNAStyle() = new{true,true}()
end

export AdaptiveStrategy
abstract type AdaptiveStrategy end
export GreedyStyle, OptimalStyle, SimpleStyle
struct GreedyStyle <: AdaptiveStrategy end
struct OptimalStyle <: AdaptiveStrategy end
struct SimpleStyle <: AdaptiveStrategy end


emptylogbook() = Array{Any,1}(undef,0)
addlogentry!(log, entry) = push!(log, entry)

# Any error criterion returns a result of the form (converged,error):
# - converged: indicates whether the approximation satisfies all convergence criteria
# - error: a criterion-specific measure of the error

# Verify the average pointwise error in a limited number of random points.
function errormeasure(::RandomPoints, platform, tolerance, f, F, args...; numrandompts = 50, verbose=false, options...)
    g = randomgrid(support(F.expansion.dictionary), numrandompts)
    z = sample(g, f)
    max_error = norm(z-F.(g), Inf)
    converged = max_error < tolerance
    verbose && @info "Errormeasure: Maximum error in $numrandompts random points did $(converged ? "" : "not ")converge "*@sprintf("%1.3e (%1.3e)",max_error,tolerance)
    converged, max_error
end

# Measure the norm of the residual of the system Ax=B.
function errormeasure(::ResidualStyle, platform, tolerance, f, F, n, A, B, C, S, L; options...)
    residual = norm(A*C-B)
    residual <= tolerance, residual
end

# Measure the norm of the residual using oversampling (as specified by the residualoversampling factor)
function errormeasure(::OversampledResidual, platform, tolerance, f, F, args...; residualoversampling=2, M = nothing, L = nothing, options...)
    # Note that we leave out any M or L parameter that the user may have passed,
    # so that we can specify an oversampling factor to the discretization routine.
    A, B = discretization(f, platform, length(F), samplingstyle=OversamplingStyle(), oversamplingfactor=residualoversampling)
    residual = norm(A*coefficients(F)-B)
    residual <= error, residual
end

# Measure the error FNA style:
# - the residual has to meet a threshold
# - the norm of the coefficients has to be smaller than a value times the estimated norm of the right hand side
function errormeasure(::FNAStyle{CREL,EREL}, platform, tolerance, f, F, n, A, B, C, S, L; optimizefase=false,
        numrandompts = 2, FNAcoef = 5.0, FNAerr = 5.0, verbose=false, Sbest=nothing, Bbest=nothing,options...) where {CREL, EREL}
    residual = norm(A*C-B)

    # Note: we pass on options, because it may contain a measure
    # Q = discrete_normalization(platform, n, L; S=S, options...)
    if optimizefase
        w = BasisFunctions.gaussweights(grid(dest(Sbest)), measure(platform))
        normF = abs(sqrt(sum(w .* Bbest.^2)))
    else
        w = BasisFunctions.gaussweights(grid(dest(S)), measure(platform))
        normF = abs(sqrt(sum(w .* B.^2)))
    end
    normErr = residual
    normCoef = norm(C)
    verbose && @info @sprintf("Errormeasure: approximate ||f|| %1.3e\n",normF)

    errConverged = EREL ? (normErr < FNAerr*normF) : (normErr < FNAerr)
    coefConverged = CREL ? (normCoef < FNAcoef*normF) : (normErr < FNAcoef)

    verbose && @info "Errormeasure: error did $(errConverged ? "" : "not ")converge: " * @sprintf("error=%1.3e, ||f||=%1.3e, c=%1.3e\n",normErr,normF,FNAerr)
    verbose && @info "Errormeasure: coefficients did $(coefConverged ? "" : "not ")converge: " * @sprintf("||c||=%1.3e, ||f||=%1.3e, c=%1.3e\n",normCoef,normF,FNAcoef)

    converged = errConverged && coefConverged
    if converged
        if optimizefase
            converged, residual
        else
            # If N is too small the norm of `f` might not be good enough. User random points for robustness.
            verbose && @info "Errormeasure: Check random points for FNA robustness. "
            converged, _ = errormeasure(RandomPoints(), platform, normF*FNAerr, f, F; numrandompts = numrandompts, verbose=verbose, options...)
            converged, residual
        end
    else
        converged, residual
    end
end


ErrorStyle(platform::Platform) = ErrorStyle(platform, DictionaryStyle(platform))
ErrorStyle(platform, ::BasisStyle) = RandomPoints()
ErrorStyle(platform, ::FrameStyle) = ResidualStyle()
ErrorStyle(platform, ::UnknownDictionaryStyle) = RandomPoints()




# Dispatch on the style of the adaptive algorithm
function approximate(fun, ap::AdaptiveApproximation;
            adaptivestyle = OptimalStyle(), criterion = ErrorStyle(ap.platform),
            verbose = false, options...)

    # Inform the user that an adaptive computation is starting
    verbose && println("\nAdaptive style: $(adaptivestyle)")
    verbose && println("Error Criterium: $criterion")

    # Do the actual computation
    F, logbook, n, tol, error, iterations, converged = adaptive_approximation(adaptivestyle, fun, ap.platform; criterion=criterion, verbose=verbose, options...)

    # Report on the final convergence
    if !converged
        @warn "Adaptive: convergence to desired tolerance not reached after $(iterations) iterations."
        if verbose
            @printf "     Tolerance: %1.3e\n" tol
            @printf "     Minimal error: %1.3e\n\n" error
        end
    elseif verbose
        @printf "     Tolerance: %1.3e\n" tol
        println("\nAdaptive: Tolerance met using $(length(F)) degrees of freedom (n=$n) in $(iterations) iterations.\n")
    end

    return F
end


function adaptive_approximation(algorithm::Symbol, args...; options...)
    if algorithm == :greedy
        style = GreedyStyle()
    elseif algorithm == :simple
        style = SimpleStyle()
    elseif algorithm == :optimal
        style = OptimalStyle()
    else
        error("Unknown adaptive algorithm.")
    end
    adaptive_approximation(style, args...; options...)
end



"""
  GreedyStyle: Greedy scheme to get decreasing coefficients.

  Let phi_i i=1..N basisfunctions.
  Take a LS of f with only phi_1, call the approximation p_1
  Take a LS of f-p_1, using only {phi_1,phi_2}, and call the approximation p_2
  Approximation f-(p_1+p_2) using {phi_k}_{k=1}^3,
  ...
  Stop the iteration when the residu of the approximation is smaller than the tolerance
  or at the maximum number of iterations.
"""
function adaptive_approximation(::GreedyStyle, f, platform;
            criterion = ErrorStyle(platform),
            p0 = param_first(platform), maxlength = 2^22, maxiterations = 100,
            threshold = 1e-12, tol=100*threshold, verbose=false, smoothing = false, options...)

    iterations = 0
    logbook = emptylogbook()

    n = p0
    F = DictFun(dictionary(platform, n))
    if smoothing
        dict = dictionary(F)
        weights = ones(length(dict))
    end

    converged = false
    while (!converged) && (length(F) <= maxlength) && (iterations <= maxiterations)
        # residual_f = x -> f(x)-F(x)
        verbose_on_first_call = verbose && (iterations==1)
        verbose_on_first_call && println("Adaptive: initial value of n is $n")
        # P, A, B, C, S, L = approximate(residual_f, platform, n; verbose=verbose_on_first_call, threshold=threshold, options...)
        # converged, error = errormeasure(criterion, platform, tol, residual_f, P, n, A, B, C, S, L; verbose=verbose, options...)
        if smoothing
            W = DiagonalOperator(dictionary(platform, n), weights)
            P, A, B, C, S, L = approximate(f, platform, n; verbose=verbose_on_first_call, threshold=threshold, D=W, options...)
        else
            P, A, B, C, S, L = approximate(f, platform, n; verbose=verbose_on_first_call, threshold=threshold, options...)
        end
        converged, error = errormeasure(criterion, platform, tol, f, P, n, A, B, C, S, L; verbose=verbose, options...)
        addlogentry!(logbook, (n, error))

        # F = F + P
        F = P

        verbose && @printf "Adaptive: N = %d, err = %1.3e, ||x|| = %1.3e\n" length(F) error norm(coefficients(F))

        n += 1
        iterations += 1
        if smoothing
            weights = vcat(weights, 1/error)
        end
    end
    return F, logbook, n, tol, error, iterations, converged
end

"SimpleStyle: the number of degrees of freedom is doubled until the tolerance is achieved."
function adaptive_approximation(::SimpleStyle, f, platform;
            criterion = ErrorStyle(platform),
            p0 = param_first(platform), maxlength = 2^22, maxiterations = 10,
            threshold = 1e-12, tol=100*threshold, verbose=false, smoothing=false, options...)

    iterations = 0
    logbook = emptylogbook()

    n = p0
    F = DictFun(dictionary(platform, n))
    if smoothing
        dict = dictionary(F)
        weights = ones(length(dict))
    end

    converged = false
    while (!converged) && (length(F) <= maxlength) && (iterations <= maxiterations)
        iterations += 1
        verbose_on_first_call = verbose && (iterations==1)
        verbose_on_first_call && println("Adaptive: initial value of n is $n")
        if smoothing
            W = DiagonalOperator(dictionary(platform, n), weights)
            F, A, B, C, S, L = approximate(f, platform, n; threshold=threshold, verbose = verbose_on_first_call, D = W, options...)
        else
            F, A, B, C, S, L = approximate(f, platform, n; threshold=threshold, verbose = verbose_on_first_call, options...)
        end
        converged, error = errormeasure(criterion, platform, tol, f, F, n, A, B, C, S, L; options...)
        addlogentry!(logbook, (n, error))

        verbose && @printf "Adaptive: N = %d, err = %1.3e, ||x|| = %1.3e\n" length(F) error norm(coefficients(F))

        nprev = n
        n = param_double(platform, n)
        if smoothing
            weights = vcat(weights, 1/error*ones(n-nprev))
        end
    end

    return F, logbook, n, tol, error, iterations, converged
end


"OptimalStyle: the number of degrees of freedom is chosen adaptively and optimally."
function adaptive_approximation(::OptimalStyle, f, platform;
        criterion = ErrorStyle(platform),
        p0 = param_first(platform), maxlength = 2^12, maxiterations = 100,
        threshold = 1e-12, tol=100*threshold, verbose=false, options...)
    iterations = 0
    logbook = emptylogbook()

    n = p0
    F = DictFun(dictionary(platform, n))

    # First let the size grow until the tolerance is reached
    A=nothing;B=nothing;C=nothing;S=nothing;L=nothing;
    previous_n = n
    next_n = n
    converged = false
    while (!converged) && (length(F) <= maxlength)
        iterations += 1
        previous_n = n
        n = next_n
        verbose_on_first_call = verbose && (iterations==1)
        verbose_on_first_call && println("Adaptive: initial value of param is $n")
        F, A, B, C, S, L = approximate(f, platform, n; threshold=threshold, verbose=verbose_on_first_call, options...)
        converged, error = errormeasure(criterion, platform, tol, f, F, n, A, B, C, S, L; verbose=verbose, options...)
        verbose && println(@sprintf("Adaptive (phase 1): N = %d, err = %1.3e, ||x|| = %1.3e ",length(F),error,norm(coefficients(F)))*" did $(converged ? "" : "not ")converge (param=$n).")

        addlogentry!(logbook, (n, error))
        next_n = param_double(platform, n)
    end

    if !converged
        # The first phase was not successful
        return F, logbook, n, tol, error, iterations, converged
    end
    verbose && println("Adaptive: Tolerance first met using $(length(F)) degrees of freedom (n=$n). Optimal n lies in [$previous_n,$n].\n")

    # Next, bisect until minimal n is found that achieves the tolerance
    lower_n = previous_n
    Fbest=F;Abest=A;Bbest=B;Cbest=C;Sbest=S;Lbest=L;
    upper_n = n

    converged2 = false
    while all(lower_n .+ 1e-4 .< upper_n)
        n = param_inbetween(platform, lower_n, upper_n)

        F, A, B, C, S, L = approximate(f, platform, n; threshold=threshold, options...)
        converged2, error = errormeasure(criterion, platform, tol, f, F, n, A, B, C, S, L;
            optimizefase=true,Sbest=Sbest,Bbest=Bbest,verbose=verbose, options...)
        addlogentry!(logbook, (n, error, lower_n, upper_n))
        if !converged2
            # At this stage, we know that lower_n does not meet the criterium, but the
            # new value of n does not either. It could be that lower_n==n, for example
            # when lower_n = 3 and upper_n = 4, and the param_inbetween method above
            # just returns 3 again. In order to avoid this, we increment n before
            # assigning to lower_n
            lower_n = param_increment(platform, n)
        else
            upper_n = n
        end
        verbose && @printf "Adaptive (phase 2): N = %d, err = %1.3e, ||x|| = %1.3e.\n" length(F) error norm(coefficients(F))
        verbose && println("Optimal param in [$(lower_n),$(upper_n)].")
        iterations += 1
    end
    if n != lower_n
        # We should take the upper_n. For simplicity, we recompute the approximation.
        n = upper_n
        F, A, B, C, S, L = approximate(f, platform, n; threshold=threshold, options...)
        converged2, error = errormeasure(criterion, platform, tol, f, F, n, A, B, C, S, L;
            optimizefase=true,Sbest=Sbest,Bbest=Bbest,verbose=verbose, options...)
    end
    verbose && println("Adaptive: Optimal n is $(n)!")

    return F, logbook, n, tol, error, iterations, converged
end

end
