
abstract type ErrorStyle end

struct ResidualStyle <: ErrorStyle end
struct OversampledResidual <: ErrorStyle end
struct RandomPoints <: ErrorStyle end
struct FNAStyle <: ErrorStyle end

abstract type AdaptiveStrategy end
struct GreedyStyle <: AdaptiveStrategy end
struct OptimalStyle <: AdaptiveStrategy end
struct SimpleStyle <: AdaptiveStrategy end


emptylogbook() = Array{Any,1}(undef,0)
addlogentry!(log, entry) = push!(log, entry)

# Any error criterion returns a result of the form (converged,error):
# - converged: indicates whether the approximation satisfies all convergence criteria
# - error: a criterion-specific measure of the error

# Verify the average pointwise error in a limited number of random points.
function errormeasure(::RandomPoints, platform, tolerance, f, F, args...; numrandompts = 50, options...)
    g = randomgrid(support(F.expansion.dictionary), numrandompts)
    z = sample(g, f)
    mean_error = sqrt(sum(abs.(z-F.(g)).^2))/numrandompts
    mean_error < tolerance, mean_error
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
function errormeasure(::FNAStyle, platform, tolerance, f, F, n, A, B, C, S, L; FNAeta = 5.0, options...)
    residual = norm(A*C-B)

    # Note: we pass on options, because it may contain a measure
    # Q = discrete_normalization(platform, n, L; S=S, options...)
    w = BasisFunctions.gaussweights(grid(dest(S)), measure(platform))
    normF = abs(sqrt(sum(w .* B.^2)))

    converged = (norm(C) < FNAeta*normF) && (residual < tolerance)
    converged, residual
end


ErrorStyle(platform::Platform) = ErrorStyle(platform, DictionaryStyle(platform))
ErrorStyle(platform, ::BasisStyle) = RandomPoints()
ErrorStyle(platform, ::FrameStyle) = ResidualStyle()
ErrorStyle(platform, ::UnknownDictionaryStyle) = RandomPoints()


# Parameter values: we define param_first, param_next and param_inbetween

# `addone(x)` adds one to x if x is a number, else it calls addone recursively on
# the elements of x
addone(x::Number) = x+1
addone(x) = map(addone, x)

# `inbetween(n1,n2)` returns an element between n1 and n2 if they are integers,
# else it calls inbetween recursively on all elements of n1 and n2.
inbetween(n1::Int, n2::Int) = (n1+n2) >> 1
inbetween(n1::T, n2::T) where {T} = map(inbetween, n1, n2)

param_first(platform::Platform) = 1
param_next(platform::Platform, n) = extension_size(dictionary(platform, n))

param_increment(platform::Platform, n) = addone(n)
param_inbetween(platform::Platform, n1, n2) = inbetween(n1, n2)


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
            threshold = 1e-12, tol=100*threshold, verbose=false, options...)

    iterations = 0
    logbook = emptylogbook()

    n = p0
    F = DictFun(dictionary(platform, n))

    converged = false
    while (!converged) && (length(F) <= maxlength) && (iterations <= maxiterations)
        residual_f = x -> f(x)-F(x)
        verbose_on_first_call = verbose && (iterations==1)
        verbose_on_first_call && println("Adaptive: initial value of n is $n")
        P, A, B, C, S = approximate(residual_f, platform, n; verbose=verbose_on_first_call, threshold=threshold, options...)
        converged, error = errormeasure(criterion, platform, tol, residual_f, P, n, A, B, C, S; options...)
        addlogentry!(logbook, (n, error))

        F = F + P

        verbose && @printf "Adaptive: N = %d, err = %1.3e, ||x|| = %1.3e\n" length(F) error norm(coefficients(F))

        n += 1
        iterations += 1
    end
    return F, logbook, n, tol, error, iterations, converged
end

"SimpleStyle: the number of degrees of freedom is doubled until the tolerance is achieved."
function adaptive_approximation(::SimpleStyle, f, platform;
            criterion = ErrorStyle(platform),
            p0 = param_first(platform), maxlength = 2^22, maxiterations = 10,
            threshold = 1e-12, tol=100*threshold, verbose=false, options...)

    iterations = 0
    logbook = emptylogbook()

    n = p0
    F = DictFun(dictionary(platform, n))

    converged = false
    while (!converged) && (length(F) <= maxlength) && (iterations <= maxiterations)
        iterations += 1
        verbose_on_first_call = verbose && (iterations==1)
        verbose_on_first_call && println("Adaptive: initial value of n is $n")
        F, A, B, C, S, L = approximate(f, platform, n; threshold=threshold, verbose = verbose_on_first_call, options...)
        converged, error = errormeasure(criterion, platform, tol, f, F, n, A, B, C, S, L; options...)
        addlogentry!(logbook, (n, error))

        verbose && @printf "Adaptive: N = %d, err = %1.3e, ||x|| = %1.3e\n" length(F) error norm(coefficients(F))

        n = param_next(platform, n)
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
    previous_n = n
    next_n = n
    converged = false
    while (!converged) && (length(F) <= maxlength)
        iterations += 1
        previous_n = n
        n = next_n
        verbose_on_first_call = verbose && (iterations==1)
        verbose_on_first_call && println("Adaptive: initial value of n is $n")
        F, A, B, C, S, L = approximate(f, platform, n; threshold=threshold, verbose=verbose_on_first_call, options...)
        converged, error = errormeasure(criterion, platform, tol, f, F, n, A, B, C, S, L; options...)
        verbose && @printf "Adaptive (phase 1): N = %d, err = %1.3e, ||x|| = %1.3e\n" length(F) error norm(coefficients(F))

        addlogentry!(logbook, (n, error))
        next_n = param_next(platform, n)
    end

    if !converged
        # The first phase was not successful
        return F, logbook, n, tol, error, iterations, converged
    end
    verbose && println("Adaptive: Tolerance first met using $(length(F)) degrees of freedom (n=$n). Optimal n lies in [$previous_n,$n].\n")

    # Next, bisect until minimal n is found that achieves the tolerance
    lower_n = previous_n
    upper_n = n

    converged2 = false
    while lower_n < upper_n
        n = param_inbetween(platform, lower_n, upper_n)

        F, A, B, C, S, L = approximate(f, platform, n; threshold=threshold, options...)
        converged2, error = errormeasure(criterion, platform, tol, f, F, n, A, B, C, S, L; options...)
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
        verbose && @printf "Adaptive (phase 2): N = %d, err = %1.3e, ||x|| = %1.3e." length(F) error norm(coefficients(F))
        verbose && println("Optimal n in [$(lower_n),$(upper_n)].")
        iterations += 1
    end
    if n != lower_n
        # We should take the upper_n. For simplicity, we recompute the approximation.
        n = upper_n
        F, A, B, C, S, L = approximate(f, platform, n; threshold=threshold, options...)
        converged2, error = errormeasure(criterion, platform, tol, f, F, n, A, B, C, S, L; options...)
    end
    verbose && println("Adaptive: Optimal n is $(n)!")

    return F, logbook, n, tol, error, iterations, converged
end
