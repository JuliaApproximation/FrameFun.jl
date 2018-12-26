
abstract type ErrorStyle end

struct RandomPoints <: ErrorStyle
end

struct OversampledResidual <: ErrorStyle
end

struct ResidualStyle <: ErrorStyle
end

function errormeasure(::RandomPoints, platform, f, F, args...; Q=50, options...)
    g = randomgrid(support(F.expansion.dictionary), Q)
    z = sqrt(sum(abs.(f.(g)-F.(g)).^2))/sum(Q)
end

function errormeasure(::ResidualStyle, platform, f, F, args...; oversamplingfactor=2, options...)
    # TODO: implement
    errormeasure(RandomPoints(), platform, f, F; options...)
end


# Generic adaptivity
function approximate(fun, ap::AdaptiveApproximation; algorithm = :optimal, options...)
    if algorithm == :greedy
        adaptive_greedy(fun, ap.platform; options...)
    elseif algorithm == :simple
        adaptive_simple(fun, ap.platform; options...)
    elseif algorithm == :optimal
        adaptive_optimal(fun, ap.platform; options...)
    else
        error("Unknown adaptive algorithm.")
    end
end

"""
  Greedy scheme to get decreasing coefficients.

  Let phi_i i=1..N basisfunctions.
  Take a LS of f with only phi_1, call the approximation p_1
  Take a LS of f-p_1, using only {phi_1,phi_2}, and call the approximation p_2
  Approximation f-(p_1+p_2) using {phi_k}_{k=1}^3,
  ...
  Stop the iteration when the residu of the approximation is smaller than the tolerance
  or at the maximum number of iterations.
"""
function adaptive_greedy(f, platform;
        criterium = ResidualStyle(),
        max_logn_coefs = 7, tol = 1e-12, verbose = false, options...)

    init_n = 4
    F = Fun(x->0, platform, init_n; options...)
    T = eltype(F)
    for n in init_n:2^max_logn_coefs
        p_i = Fun(x->(f(x)-F(x)), platform, n; coefficienttype = T, options...)
        F = F + p_i
        R = errormeasure(criterium, platform, f, F)
        verbose && println("Adaptive: using $n degrees of freedom, residual $R")
        if R < tol
            verbose && println("Adaptive: stopped with residual $R")
            return F
        end
    end
    F
end

"""
  Create approximation to function with a given dictionary in a domain.

  The number of points is chosen adaptively, and guarantees the tolerance, but may not be optimal.
"""
function adaptive_simple(f::Function, platform;
        no_checkpoints = 50, max_logn_coefs=8, tol=1e-12, abscoef=nothing,
        verbose=false, options...)

    n = 4
    dict = Dictionary(platform, n)
    domain = support(dict)
    T = codomaintype(dict)
    N = dimension(domain)
    F = nothing
    rgrid = randomgrid(domain, no_checkpoints)
    error = Inf
    random_f = sample(rgrid, f, eltype(f(rgrid[1]...)))
    random_F = zeros(T,no_checkpoints)

    while length(dict) <= 1<<(max_logn_coefs)
        n = extension_size(dict)
        dict = Dictionary(platform, n)
        F = Fun(f, platform, n; options...)
        random_F = F(rgrid)
        error = maximum(abs.(random_F-random_f))
        verbose && (@printf "Adaptive: error with %d coefficients is %1.3e\n" (length(dict)) error)
        (verbose & (abscoef != nothing))  && (@printf "Norm with %d coefficients is %1.3e (%1.3e)\n" (length(dict)) norm(coefficients(F)) abscoef*_trapnorm(f, dictionary(F), 1))
        if (error < tol) && ((abscoef == nothing) || abs_coefficient_test(f, F, abscoef, 1))
            return F
        end
    end
    @warn("Adaptive: maximum number of coefficients exceeded, error is $(error)")
    F
end

"""
  Create approximation to function with a dictionary in a domain.

  The number of points is chosen adaptively and optimally.
"""
function adaptive_optimal(f::Function, platform;
        criterium = ResidualStyle(),
        no_checkpoints=50, max_logn_coefs=8, threshold = 1e-10, tol=100*threshold, verbose=false, adaptive_verbose = verbose, return_log=false, randomtest=false, abscoef=nothing, relcoef=nothing, options...)

    coefficienttest = (nothing!=abscoef) | (nothing!=relcoef)

    n = 8
    dict = Dictionary(platform, n)
    domain = support(dict)
    ELT = codomaintype(dict)
    N = dimension(dict)
    F = nothing
    rgrid = randomgrid(domain, no_checkpoints)
    Nmax = NaN; Nmin = 1;
    return_log && (log = zeros(T,0,4))
    its = 0
    error = -1
    while length(dict) <= 2^max_logn_coefs && its < 100
        # Find new approximation
        # F, A, B, C, D, S = approximate(f, dict, domain; threshold=threshold, options...)
        F = Fun(f, platform, n; threshold=threshold, options...)

        # Using residual
        error = errormeasure(criterium, platform, f, F; options...)
        # println(n, ": ", Nmin, " ", Nmax, " ", error)
        return_log && (log = [log; n Nmin Nmax error])
        adaptive_verbose  && (@printf "Adaptive: error with %d coefficients is %1.3e (%1.3e)\n" (length(dict)) error tol)

        errorsucces = (error < tol)

        randomsuccess = randomtest ?
            random_test(f, F, tol, 3) : true
        (adaptive_verbose&randomtest)  && (randomsuccess ?
            (@printf "Adaptive: random test in 3 points succeeded\n") : (@printf "Random test in 3 points failed\n"))

        coefsuccess = coefficienttest ?
            coefficient_test(f, F; tolerance=tol, abscoef=abscoef, relcoef=relcoef, options...) : true
        (adaptive_verbose&coefficienttest)  && (coefsuccess ?
            (@printf "Adaptive: coefficient test succeeded %1.3e\n" norm(coefficients(F))) :
            (@printf "Adaptive: coefficient test failed %1.3e\n" norm(coefficients(F))))

        success = broadcast(&, errorsucces, coefsuccess, randomsuccess)
        success ? Nmax=n : Nmin=n

        # Stop condition
        # If the bounds Nmin and Nmax are determined and if they are close
        if are_close(Nmin, Nmax)
            dict = Dictionary(platform, Nmax)
            F = Fun(f, platform, Nmax; threshold=threshold, options...)

            return_log && (return F, log)
            return F
        end

        # Decide on new n
        if isnan(Nmax)
            # tolerance is not yet met, increase n
            n = extension_size(dict)
        else
            # tolerance is reached, find optimal n by subdivision
            if N==1
                # in one dimension Nmin and Nmax are floats
                n = round(Int,(Nmin + Nmax)/2)
            else
                # in multiple dimensions, Nmin and Nmax are tupples
                n = (round.(Int,(collect(Nmin) + collect(Nmax))/2)...,)
            end
        end

        # Use new n
        dict = Dictionary(platform, n)
        its = its + 1
    end
    @warn("Adaptive: maximum number of coefficients exceeded, error is $(error)")
    F
end
