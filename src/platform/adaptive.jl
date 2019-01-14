
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

# Any error criterium returns a result of the form (converged,error):
# - converged: indicates whether the approximation satisfies all convergence criteria
# - error: a criterium-specific measure of the error

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
function errormeasure(::FNAStyle, platform, tolerance, f, F, n, A, B, C, S, L; eta = 5.0, options...)
    residual = norm(A*C-B)

    Q = discrete_normalization(platform, n, L; S=S)
    normF = abs(sqrt(sum(Q * B.^2)))

    converged = (norm(C) < eta*normF) && (residual < tolerance)
    converged, residual
end


nextsize(platform::Platform, n; dict = dictionary(platform, n)) = extension_size(dict)

ErrorStyle(platform::Platform) = ErrorStyle(platform, DictionaryStyle(platform))
ErrorStyle(platform, ::BasisStyle) = RandomPoints()
ErrorStyle(platform, ::FrameStyle) = ResidualStyle()
ErrorStyle(platform, ::UnknownDictionaryStyle) = RandomPoints()

# Initial values

first_parameters(p::Platform) = (8,8)
first_sizeparameter(p::Platform) = first_parameters(p)[1]
first_samplingparameter(p::Platform) = first_parameters(p)[2]



# Dispatch on the style of the adaptive algorithm
function approximate(fun, ap::AdaptiveApproximation;
            adaptivestyle = OptimalStyle(), criterium = ErrorStyle(ap.platform),
            verbose = false, options...)

    # Inform the user that an adaptive computation is starting
    verbose && println("\nAdaptive style: $(adaptivestyle)")
    verbose && println("Error Criterium: $criterium")

    # Do the actual computation
    F, logbook, n, tol, error, iterations, converged = adaptive_approximation(adaptivestyle, fun, ap.platform; criterium=criterium, verbose=verbose, options...)

    # Report on the final convergence
    if !converged
        @warn "Adaptive: convergence to desired tolerance not reached after $(iterations) iterations."
        if verbose
            @printf "     Tolerance: %1.3e\n" tol
            @printf "     Minimal error: %1.3e\n\n" error
        end
    elseif verbose
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
            criterium = ErrorStyle(platform),
            initial_n = 1, maxlength = 2^22, maxiterations = 100,
            threshold = 1e-12, tol=100*threshold, verbose=false, options...)

    local error
    iterations = 0
    logbook = emptylogbook()

    n = initial_n
    F = DictFun(dictionary(platform, n))

    converged = false
    while (!converged) && (length(F) <= maxlength) && (iterations <= maxiterations)
        residual_f = x -> f(x)-F(x)
        verbose_on_first_call = verbose && (iterations==1)
        verbose_on_first_call && println("Adaptive: initial value of n is $n")
        P, A, B, C, S = approximate(residual_f, platform, n; verbose=verbose_on_first_call, threshold=threshold, options...)
        converged, error = errormeasure(criterium, platform, tol, residual_f, P, n, A, B, C, S; options...)
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
            criterium = ErrorStyle(platform),
            initial_n = 8, maxlength = 2^22, maxiterations = 10,
            threshold = 1e-12, tol=100*threshold, verbose=false, options...)

    local error
    iterations = 0
    logbook = emptylogbook()

    n = initial_n
    F = DictFun(dictionary(platform, n))

    converged = false
    while (!converged) && (length(F) <= maxlength) && (iterations <= maxiterations)
        iterations += 1
        verbose_on_first_call = verbose && (iterations==1)
        verbose_on_first_call && println("Adaptive: initial value of n is $n")
        F, A, B, C, S, L = approximate(f, platform, n; threshold=threshold, verbose = verbose_on_first_call, options...)
        converged, error = errormeasure(criterium, platform, tol, f, F, n, A, B, C, S, L; options...)
        addlogentry!(logbook, (n, error))

        verbose && @printf "Adaptive: N = %d, err = %1.3e, ||x|| = %1.3e\n" length(F) error norm(coefficients(F))

        n = nextsize(platform, n)
    end

    return F, logbook, n, tol, error, iterations, converged
end


"OptimalStyle: the number of degrees of freedom is chosen adaptively and optimally."
function adaptive_approximation(::OptimalStyle, f, platform;
        criterium = ErrorStyle(platform),
        initial_n = 8, maxlength = 2^12, maxiterations = 100,
        threshold = 1e-12, tol=100*threshold, verbose=false, options...)

    local error
    iterations = 0
    logbook = emptylogbook()

    n = initial_n
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
        converged, error = errormeasure(criterium, platform, tol, f, F, n, A, B, C, S, L; options...)
        verbose && @printf "Adaptive (phase 1): N = %d, err = %1.3e, ||x|| = %1.3e\n" length(F) error norm(coefficients(F))

        addlogentry!(logbook, (n, error))
        next_n = nextsize(platform, n)
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
        n = (upper_n+lower_n) >> 1

        F, A, B, C, S, L = approximate(f, platform, n; threshold=threshold, options...)
        converged2, error = errormeasure(criterium, platform, tol, f, F, n, A, B, C, S, L; options...)
        addlogentry!(logbook, (n, error, lower_n, upper_n))
        if !converged2
            lower_n = n+1
        else
            upper_n = n
        end
        verbose && @printf "Adaptive (phase 2): N = %d, err = %1.3e, ||x|| = %1.3e. Optimal n in [%d,%d].\n" length(F) error norm(coefficients(F)) lower_n upper_n
        iterations += 1
    end
    if n != lower_n
        # We should take the upper_n. For simplicity, we recompute the approximation.
        n = upper_n
        F, A, B, C, S, L = approximate(f, platform, n; threshold=threshold, options...)
        converged2, error = errormeasure(criterium, platform, tol, f, F, n, A, B, C, S, L; options...)
    end
    verbose && println("Adaptive: Optimal n is $(n)!")

    return F, logbook, n, tol, error, iterations, converged
end
