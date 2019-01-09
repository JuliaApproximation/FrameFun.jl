
abstract type ErrorStyle end

struct RandomPoints <: ErrorStyle end
struct OversampledResidual <: ErrorStyle end
struct ResidualStyle <: ErrorStyle end

abstract type AdaptiveStrategy end
struct GreedyStyle <: AdaptiveStrategy end
struct OptimalStyle <: AdaptiveStrategy end
struct SimpleStyle <: AdaptiveStrategy end


emptylogbook() = Array{Any,1}(undef,0)
addlogentry!(log, entry) = push!(log, entry)


function errormeasure(::RandomPoints, platform, f, F, args...; numrandompts = 50, options...)
    g = randomgrid(support(F.expansion.dictionary), numrandompts)
    z = sqrt(sum(abs.(f.(g)-F.(g)).^2))/numrandompts
end

errormeasure(::ResidualStyle, platform, f, F, A, B, C, S; oversamplingfactor=2, options...) = norm(A*C-B)


nextsize(platform::Platform, n; dict = dictionary(platform, n)) = extension_size(dict)

ErrorStyle(platform::Platform) = ErrorStyle(platform, DictionaryStyle(platform))
ErrorStyle(platform, ::BasisStyle) = RandomPoints()
ErrorStyle(platform, ::FrameStyle) = ResidualStyle()
ErrorStyle(platform, ::UnknownDictionaryStyle) = RandomPoints()


# Dispatch on the style of the adaptive algorithm
function approximate(fun, ap::AdaptiveApproximation; adaptivestyle = OptimalStyle(), verbose = false, options...)
    verbose && println("\nAdaptive style: $(adaptivestyle)")

    F, logbook, n, tol, error, iterations = adaptive_approximation(adaptivestyle, fun, ap.platform; verbose=verbose, options...)

    if error > tol
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

    n = initial_n
    F = DictFun(dictionary(platform, n))
    error = typemax(domaintype(dictionary(F)))

    iterations = 0
    logbook = emptylogbook()

    while (error > tol) && (length(F) <= maxlength) && (iterations <= maxiterations)
        residual_f = x -> f(x)-F(x)
        P, A, B, C, S = approximate(residual_f, platform, n; threshold=threshold, options...)
        error = errormeasure(criterium, platform, residual_f, P, A, B, C, S; options...)
        addlogentry!(logbook, (n, error))

        F = F + P

        verbose && @printf "Adaptive: error with %d coefficients is %1.3e (tolerance: %1.3e)\n" length(F) error tol

        n += 1
        iterations += 1
    end
    return F, logbook, n, tol, error, iterations
end

"SimpleStyle: the number of degrees of freedom is doubled until the tolerance is achieved."
function adaptive_approximation(::SimpleStyle, f, platform;
            criterium = ErrorStyle(platform),
            initial_n = 1, maxlength = 2^22, maxiterations = 10,
            threshold = 1e-12, tol=100*threshold, verbose=false, options...)

    logbook = emptylogbook()
    iterations = 0

    n = initial_n
    F = DictFun(dictionary(platform, n))
    error = typemax(domaintype(dictionary(F)))

    while (error > tol) && (length(F) <= maxlength) && (iterations <= maxiterations)
        F, A, B, C, S = approximate(f, platform, n; threshold=threshold, options...)
        error = errormeasure(criterium, platform, f, F, A, B, C, S; options...)
        addlogentry!(logbook, (n, error))

        verbose && @printf "Adaptive: error with %d coefficients is %1.3e (tolerance: %1.3e)\n" length(F) error tol

        n = nextsize(platform, n)
        iterations += 1
    end

    return F, logbook, n, tol, error, iterations
end


"OptimalStyle: the number of degrees of freedom is chosen adaptively and optimally."
function adaptive_approximation(::OptimalStyle, f, platform;
        criterium = ErrorStyle(platform),
        initial_n = 1, maxlength = 2^12, maxiterations = 100,
        threshold = 1e-12, tol=100*threshold, verbose=false, options...)

    logbook = emptylogbook()
    iterations = 0

    n = initial_n
    F = DictFun(dictionary(platform, n))
    error = typemax(domaintype(dictionary(F)))

    # First let the size grow until the tolerance is reached
    previous_n = n
    next_n = n
    while (error > tol) && (length(F) <= maxlength)
        previous_n = n
        n = next_n
        F, A, B, C, S = approximate(f, platform, n; threshold=threshold, options...)
        error = errormeasure(criterium, platform, f, F, A, B, C, S; options...)
        verbose && @printf "Adaptive (phase 1): error with %d coefficients is %1.3e (tolerance: %1.3e)\n" length(F) error tol

        addlogentry!(logbook, (n, error))
        next_n = nextsize(platform, n)
        iterations += 1
    end

    if error > tol
        # The first phase was successful
        return F, logbook, n, tol, error, iterations
    end
    verbose && println("Adaptive: Tolerance first met using $(length(F)) degrees of freedom (n=$n). Optimal n lies in [$previous_n,$n].\n")

    # Next, bisect until minimal n is found that achieves the tolerance
    lower_n = previous_n
    upper_n = n

    while lower_n < upper_n
        n = (upper_n+lower_n) >> 1

        F, A, B, C, S = approximate(f, platform, n; threshold=threshold, options...)
        error = errormeasure(criterium, platform, f, F, A, B, C, S; options...)
        addlogentry!(logbook, (n, error, lower_n, upper_n))
        if error > tol
            lower_n = n+1
        else
            upper_n = n
        end
        if verbose
            @printf "Adaptive (phase 2): n = %d, error: %1.3e. Optimal n in [%d,%d].\n" n error lower_n upper_n
        end
        iterations += 1
    end
    if n != lower_n
        # We should take the upper_n. For simplicity, we recompute the approximation.
        n = upper_n
        F, A, B, C, S = approximate(f, platform, n; threshold=threshold, options...)
        error = errormeasure(criterium, platform, f, F, A, B, C, S; options...)
    end
    (error <= tol) && verbose && println("Adaptive: Optimal n is $(n)!")

    return F, logbook, n, tol, error, iterations
end
