module Adaptivity
using Printf, FrameFun.Platforms, BasisFunctions, FrameFun.ApproximationProblems

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
struct OptimalStyleFirstFase <: AdaptiveStrategy end


emptylogbook() = Array{Any,1}(undef,0)
addlogentry!(log, entry) = push!(log, entry)

# Any error criterion returns a result of the form (converged,error):
# - converged: indicates whether the approximation satisfies all convergence criteria
# - error: a criterion-specific measure of the error

# Verify the average pointwise error in a limited number of random points.
function errormeasure(::RandomPoints, platform, tolerance, f, F, args...; numrandompts = 50, verbose=false, options...)
    g = randomgrid(support(F), numrandompts)
    z = sample(g, f)
    max_error = norm(z-F(g), Inf)
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


function adaptive_approximate_step!(logbook, f, platform, n, δ; criterion = ErrorStyle(platform), options...)
    F, A, B, C, S, L = approximate(f, platform, n; options...)
    converged, error = errormeasure(criterion, platform, δ, f, F, n, A, B, C, S, L; options...)
    addlogentry!(logbook, (n, error))
    F, A, B, C, S, L, converged, error
end

function initial_weight(platform, nprev, n, error)
    dict = dictionary(platform, n)
    ScalingOperator(dict, error)
end

function restrict_weight(platform, W, nbest, n)
    dictbest = dictionary(platform, nbest)
    dict = dictionary(platform, n)
    DiagonalOperator(dict, restriction(dictbest, dict)*diag(W))
end

function next_weight(platform, W, nprev, n, error)
    dict_small = dictionary(platform, nprev)
    dict_large = dictionary(platform, n)
    E = extension(dict_small, dict_large)
    D = E*diag(W)
    for i in eachindex(D)
        # might not be robust
        if D[i] == 0
            D[i] = error
        end
    end
    DiagonalOperator(dict_large, D)
end

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


# "SimpleStyle: the number of degrees of freedom is doubled until the tolerance is achieved."
# "GreedyStyle: Greedy scheme to get decreasing coefficients."
for (STYLE,nextparam,ret) in zip((:GreedyStyle,:SimpleStyle,:OptimalStyleFirstFase),(:param_increment,:param_double,:param_double),
        (Meta.parse("F, logbook, n, δ, error, iterations, converged"),
        Meta.parse("F, logbook, n, δ, error, iterations, converged"),
        Meta.parse("F, A, B, C, S, L, logbook, n, nprev, δ, error, iterations, converged, (weightedAZ ? (W) : ())...")) )
    @eval function adaptive_approximation(::$STYLE, f, platform;
                p0 = param_first(platform), maxlength = 2^22, maxiterations = 10,
                threshold = 1e-12, δ=100*threshold, verbose=false, weightedAZ=false, options...)

        iterations = 0
        logbook = emptylogbook()
        n = p0
        nprev = p0
        # First approximation
        verbose && println("Adaptive: initial value of n is $n")
        iterations += 1
        F, A, B, C, S, L, converged, error = adaptive_approximate_step!(logbook, f, platform, n, δ; threshold=threshold, verbose=verbose, options...)
        verbose && @printf "Adaptive: N = %d, err = %1.3e, ||x|| = %1.3e\n" length(F) error norm(coefficients(F))
        converged && (return $ret)

        # Preparation for next step
        nprev = n
        n = $nextparam(platform, n)
        if weightedAZ
            W = initial_weight(platform, nprev, n, error)
        end

        # Next Steps
        while (length(F) <= maxlength) && (iterations <= maxiterations)
            iterations += 1

            # If weightedAZ add the weight
            if weightedAZ
                F, A, B, C, S, L, converged, error = adaptive_approximate_step!(logbook, f, platform, n, δ;
                    threshold=threshold, verbose=false,
                    weightedAZ = true, AZ_Cweight = W,
                    options...)
            else
                F, A, B, C, S, L, converged, error = adaptive_approximate_step!(logbook, f, platform, n, δ;
                    threshold=threshold, verbose=false,
                    options...)
            end
            verbose && @printf "Adaptive: N = %d, err = %1.3e, ||x|| = %1.3e\n" length(F) error norm(coefficients(F))
            converged && (return $ret)

            # Preparation for next step
            nprev = n
            n = $nextparam(platform, n)
            if weightedAZ
                W = next_weight(platform, W, nprev, n, error)
            end
        end

        return $ret
    end
end

# "OptimalStyle: the number of degrees of freedom is chosen adaptively and optimally."
function adaptive_approximation(::OptimalStyle, f, platform;
        maxlength = 2^12, maxiterations = 100,
        threshold = 1e-12, δ=100*threshold, verbose=false, weightedAZ=false, stoptolerance=δ, options...)

    # First let the size grow until the tolerance is reached
    R = adaptive_approximation(OptimalStyleFirstFase(), f, platform;
            weightedAZ = weightedAZ, maxlength = maxlength, maxiterations = maxiterations,
            threshold = threshold, δ=δ, verbose=verbose, options...)

    if weightedAZ
        Fbest, Abest, Bbest, Cbest, Sbest, Lbest, logbook, upper_n, lower_n, δ, error, iterations, converged, W = R
    else
        Fbest, Abest, Bbest, Cbest, Sbest, Lbest, logbook, upper_n, lower_n, δ, error, iterations, converged = R
    end
    n = nbest = upper_n
    F, A, B, C, S, L = Fbest, Abest, Bbest, Cbest, Sbest, Lbest

    # The first phase was not successful
    if !converged
        verbose && println("Adaptive: The first fase was not successful. Optimal n lies beyond $nbest.\n")
        return Fbest, logbook, upper_n, δ, error, iterations, converged
    end

    # The first phase was successful
    verbose && println("Adaptive: Tolerance first met using $(length(Fbest)) degrees of freedom (n=$nbest). Optimal n lies in [$lower_n,$upper_n].\n")


    # Next, bisect until minimal n is found that achieves the tolerance
    while hasparam_inbetween(platform, lower_n, upper_n, stoptolerance)
        n = param_inbetween(platform, lower_n, upper_n)

        if weightedAZ
            F, A, B, C, S, L, converged, error = adaptive_approximate_step!(logbook, f, platform, n, δ;
                threshold=threshold, verbose=false,
                weightedAZ = true, AZ_Cweight = restrict_weight(platform, W, nbest, n),
                optimizefase=true,Sbest=Sbest,Bbest=Bbest,
                options...)
        else
            F, A, B, C, S, L, converged, error = adaptive_approximate_step!(logbook, f, platform, n, δ;
                threshold=threshold, verbose=false,
                optimizefase=true,Sbest=Sbest,Bbest=Bbest,
                options...)
        end

        if !converged
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
        if weightedAZ
            F, A, B, C, S, L, converged, error = adaptive_approximate_step!(logbook, f, platform, n, δ;
                threshold=threshold, verbose=false,
                weightedAZ = true, AZ_Cweight = restrict_weight(platform, W, nbest, n),
                optimizefase=true,Sbest=Sbest,Bbest=Bbest,
                options...)
        else
            F, A, B, C, S, L, converged, error = adaptive_approximate_step!(logbook, f, platform, n, δ;
                threshold=threshold, verbose=false,
                optimizefase=true,Sbest=Sbest,Bbest=Bbest,
                options...)
        end

    end
    verbose && println("Adaptive: Optimal n is $(n)!")

    return F, logbook, n, δ, error, iterations, converged
end


end
