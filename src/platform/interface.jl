##########################
# Approximation problems
##########################

"""
Approximation problem types group the arguments that a user supplies
to the `Fun` constructor. This may be a dictionary or a platform, for example.

A concrete approximaton problem implements an interface. For example, one can
ask for its sampling operator. However, all functionality is delegated to the
underlying dictionary or platform. Thus, an approximation problem is nothing but
an empty intermediate layer that passes on the questions to the right object.

In turn, by default a platform delegates all queries to the dictionary it
represents. Hence, all defaults in the interface are really specified by the
dictionary. The platform can intercept and change anything, for example in order
to specify a sampling operator that is different from the default of the dictionary.

Other routines besides `Fun` can group their arguments into an approximation
problem too. This ensures that these routines automatically implement the exact
same interface as the `Fun` constructor. Thus, it is easy to write a
`samplingoperator(...)` routine that accepts the exact same arguments as `Fun`,
and that simply returns the sampling operator that the `Fun` constructor would
use to solve the approximation problem.
"""
abstract type ApproximationProblem end

"""
A `DictionaryApproximation` stores a concrete dictionary that was supplied to
the Fun constructor. All queries are directed to the dictionary and hence they
are all defaults.

The user may override the defaults by explicitly specifying optional arguments
to the `Fun` constructor, for example `samplingstyle=InterpolationStyle()`. Such
optional arguments overrule the defaults of the dictionary.
"""
mutable struct DictionaryApproximation <: ApproximationProblem
    dict    ::  Dictionary
    samplingparam

    DictionaryApproximation(dict::Dictionary, samplingparam = nothing) =
        new(dict, samplingparam)
end

dictionary(ap::DictionaryApproximation) = ap.dict

coefficienttype(ap::DictionaryApproximation) = coefficienttype(ap.dict)

elements(ap::DictionaryApproximation) = map(approximationproblem, elements(ap.dict))


"""
A `PlatformApproximation` stores a platform and a parameter value. This
corresponds to a specific dictionary (which is also stored). All queries are
delegated to the platform, and thus the answers may differ from the defaults.

As with dictionary approximations, the settings can be overruled by explicitly
specifying optional arguments to `Fun`.
"""
mutable struct PlatformApproximation <: ApproximationProblem
    platform    ::  Platform
    dict        ::  Dictionary
    param
    samplingparam

    PlatformApproximation(platform, param, samplingparam = nothing) =
        new(platform, dictionary(platform, param), param, samplingparam)
end

dictionary(ap::PlatformApproximation) = ap.dict
platform(ap::PlatformApproximation) = ap.platform
parameter(ap::PlatformApproximation) = ap.param

coefficienttype(ap::PlatformApproximation) = coefficienttype(ap.dict)



"""
An `AdaptiveApproximation` only stores a platform, with which to compute adaptive
approximations.
"""
struct AdaptiveApproximation <: ApproximationProblem
    platform    ::  Platform
end

dictionary(ap::AdaptiveApproximation, param) = ap.platform[param]


## The construction of an approximation problem

# 1. From a dictionary:
approximationproblem(dict::Dictionary) = approximationproblem(coefficienttype(dict), dict)
approximationproblem(dict::Dictionary, domain::Domain) =
    approximationproblem(promote_type(coefficienttype(dict),eltype(domain)), dict, domain)
# This two-argument routine allows to specify a coefficient type, in which case the
# dictionary will be promoted if necessary. This allows the `Fun` constructor to ensure
# that the dictionary can handle complex coefficients, for example.
approximationproblem(::Type{T}, dict::Dictionary) where {T} =
    DictionaryApproximation(promote_coefficienttype(dict, T))

# If a dictionary and a domain is specified, we make an extension frame.
function approximationproblem(::Type{T}, dict::Dictionary, domain::Domain) where {T}
    if domain == support(dict)
        approximationproblem(T, dict)
    else
        approximationproblem(T, extensionframe(domain, promote_coefficienttype(dict, T)))
    end
end


# 2. An approximation problem from a platform and a concrete parameter value
approximationproblem(platform, param) = PlatformApproximation(platform, param)
approximationproblem(platform, param, L) = PlatformApproximation(platform, param, L)

# 3. An adaptive approximation problem from a platform
approximationproblem(platform) = AdaptiveApproximation(platform)

# From the adaptive approximation problem we can create a concrete platform
# approximation by specifying a parameter value.
approximationproblem(ap::AdaptiveApproximation, param) = approximationproblem(ap.platform, param)
approximationproblem(ap::AdaptiveApproximation, param, L) = approximationproblem(ap.platform, param, L)



################################################
# Routines that accept the Fun interface
################################################

# Below is an exhaustive list of functions that implement the Fun interface.

# The sampling and dual sampling operator
for op in (:samplingoperator, :dualsamplingoperator, :samplingparameter, :sampling_grid)
    @eval $op(dict::Dictionary, args...; options...) = $op(approximationproblem(dict, args...); options...)
    @eval $op(platform::Platform, args...; options...) = $op(approximationproblem(platform, args...); options...)
end

# The discretization and dualdiscretization
for op in (:discretization, :dualdiscretization)
    @eval $op(dict::Dictionary, args...; options...) = $op(approximationproblem(dict, args...); options...)
    @eval $op(platform::Platform, args...; options...) = $op(approximationproblem(platform, args...); options...)
end

# The discretization routine can also take f as an argument
discretization(f, dict::Dictionary, args...; options...) = discretization(f, approximationproblem(dict, args...); options...)
discretization(f, platform::Platform, args...; options...) = discretization(f, approximationproblem(platform, args...); options...)

# The solvers and related operators
for op in (:solver, :AZ_A, :AZ_Z, :AZ_Zt, :plungeoperator, :plungematrix, :plungerank, :smoothingoperator)
    @eval $op(dict::Dictionary, args...; options...) = $op(approximationproblem(dict, args...); options...)
    @eval $op(platform::Platform, args...; options...) = $op(approximationproblem(platform, args...); options...)
end



########################
# Platform delegation
########################

# We delegate some functionality from the platform to the dictionary corresponding
# to a value of the platform parameter n.
# Platforms may choose to override these to enable different behaviour.

interpolation_grid(platform::Platform, n; dict = dictionary(platform, n), options...) =
    interpolation_grid(dict)

oversampling_grid(platform::Platform, n, L; dict = dictionary(platform, n), options...) =
    oversampling_grid(dict, L)

dualdictionary(platform::Platform, n; dict = dictionary(platform, n)) =
    dualdictionary(dict)

discrete_normalization(platform::Platform, n, L; S = samplingoperator(platform, n, L), options...) =
    quadraturenormalization(S, measure(platform))

space(platform::Platform) = space(measure(platform))


########################
# Interface delegation
########################

# Here, we direct queries to the dictionary or to the platform.
# Also, queries to the platform are then redirected to the dictionary. Some
# of the redirection is implement in platform.jl and not here.

## SamplingStyle and SolverStyle

SamplingStyle(ap::DictionaryApproximation) = SamplingStyle(dictionary(ap))
SolverStyle(ap::DictionaryApproximation, dstyle::SamplingStyle) = SolverStyle(dictionary(ap), dstyle)

SamplingStyle(ap::PlatformApproximation) = SamplingStyle(ap.platform)
SolverStyle(ap::PlatformApproximation, dstyle::SamplingStyle) = SolverStyle(ap.platform, dstyle)


## The sampling parameter

samplingparameter(ap::ApproximationProblem; samplingstyle = SamplingStyle(ap), options...) =
    samplingparameter(samplingstyle, ap; options...)

function samplingparameter(samplingstyle::SamplingStyle, ap::ApproximationProblem; options...)
    # Has it been computed before or was it supplied by the user?
    if ap.samplingparam != nothing
        ap.samplingparam
    else
        # It wasn't. We deduce its value from the options given and store the outcome.
        L = deduce_samplingparameter(samplingstyle, ap; options...)
        ap.samplingparam = L
    end
end

function deduce_samplingparameter(::InterpolationStyle, ap;
            verbose = false, oversamplingfactor = 2, options...)
    N = length(dictionary(ap))
    L = match_sampling_parameter(ap, N)
    L
end

function deduce_samplingparameter(::OversamplingStyle, ap;
            verbose = false, oversamplingfactor = 2, options...)
    if haskey(options, :L)
        # The user specified L as an option
        L = options[:L]
        verbose && println("Sampling parameter: using L = $L")
        return L
    end
    # In the absence of L, we deduce M and then find the best matching L
    # M is either supplied, or we compute it based on the (default) oversamplingfactor
    M = haskey(options, :M) ? options[:M] : round(Int, oversamplingfactor * length(dictionary(ap)))
    L = match_sampling_parameter(ap, M)
    verbose && println("Sampling parameter: best match for M = $M is L = $L")
    L
end

# TODO: implement this one better (more general)
deduce_samplingparameter(::GramStyle, ap; options...) = length(dictionary(ap))

deduce_samplingparameter(::RectangularGramStyle, ap; projectiondict, options...) = length(projectiondict)

deduce_samplingparameter(::GridStyle, ap; options...) = nothing

samplingparameter(samplingstyle::ProductSamplingStyle, ap::ApproximationProblem; options...) =
    map((x,style)->samplingparameter(x; samplingstyle=style, options...), elements(ap), samplingstyle.styles)


## Sampling operator

# We dispatch on the sampling style
samplingoperator(ap::ApproximationProblem; samplingstyle = SamplingStyle(ap), options...) =
    samplingoperator(samplingstyle, ap; options...)

dualsamplingoperator(ap::ApproximationProblem; samplingstyle = SamplingStyle(ap), options...) =
    dualsamplingoperator(samplingstyle, ap; options...)


function samplingoperator(samplingstyle::DiscreteStyle, ap::ApproximationProblem;
            T = coefficienttype(ap), normalizedsampling = false, options...)
    grid = sampling_grid(samplingstyle, ap; options...)
    if normalizedsampling
        D = sampling_normalization(GridBasis{T}(grid), measure(ap); T=T, options...)
        D * GridSampling(grid, T)
    else
        GridSampling(grid, T)
    end
end

samplingoperator(::GenericSamplingStyle, ap::PlatformApproximation; options...) =
    genericsamplingoperator(ap.platform, ap.param; dict = ap.dict, options...)

measure(ap::DictionaryApproximation) = measure(dictionary(ap))
measure(ap::PlatformApproximation) = measure(platform(ap))

samplingoperator(samplingstyle::GramStyle, ap::ApproximationProblem;
            measure = measure(ap), options...) =
    ProjectionSampling(dictionary(ap), measure)

samplingoperator(samplingstyle::RectangularGramStyle, ap::ApproximationProblem;
            projectiondict, measure = measure(ap), options...) =
    ProjectionSampling(projectiondict, measure)

samplingoperator(samplingstyle::ProductSamplingStyle, ap::ApproximationProblem; options...) =
    tensorproduct( map( (x,style) -> samplingoperator(x; samplingstyle=style, options...),
                        elements(ap), samplingstyle.styles)...
    )

# First, if necessary, we recompute the samplingoperator using the given options.
function dualsamplingoperator(samplingstyle::DiscreteStyle, ap::ApproximationProblem; options...)
    local S
    if haskey(options, :S)
        S = options[:S]
    else
        S = samplingoperator(samplingstyle, ap; options...)
    end
    dualsamplingoperator(samplingstyle, ap, S; options...)
end

function dualsamplingoperator(samplingstyle::DiscreteStyle, ap::ApproximationProblem, S;
            normalizedsampling = false, options...)
    if normalizedsampling
        S
    else
        O = discrete_normalization(ap; S=S)
        O * S
    end
end

dualsamplingoperator(samplingstyle::ProductSamplingStyle, ap::ApproximationProblem, S;
            options...) =
    tensorproduct(
        map( (x,Sel,style) -> dualsamplingoperator(x; S=Sel, samplingstyle=style, options...),
             elements(ap), productelements(S), samplingstyle.styles)...
    )

## Discrete sampling grid

sampling_grid(ap::ApproximationProblem; samplingstyle = SamplingStyle(ap), options...) =
    sampling_grid(samplingstyle, ap; options...)

# - interpolation: we invoke interpolation_grid on the dictionary or platform
sampling_grid(::InterpolationStyle, ap::DictionaryApproximation; options...) =
    interpolation_grid(dictionary(ap))
sampling_grid(::InterpolationStyle, ap::PlatformApproximation; options...) =
    interpolation_grid(ap.platform, ap.param; dict = ap.dict, options...)

# - generic grid: we invoke platform_grid on the platform
sampling_grid(::GridStyle, ap::DictionaryApproximation; grid, options...) = grid
sampling_grid(::GridStyle, ap::PlatformApproximation; options...) =
    platform_grid(ap.platform, ap.param; dict = ap.dict, options...)
platform_grid(platform::Platform, param; grid, options...) = grid

# - oversampling: we invoke oversampling_grid on the dictionary or platform
oversampling_grid(ap::DictionaryApproximation, L; options...) =
    oversampling_grid(dictionary(ap), L)
oversampling_grid(ap::PlatformApproximation, L; options...) =
    oversampling_grid(ap.platform, ap.param, L; dict = dictionary(ap), options...)

function sampling_grid(sstyle::OversamplingStyle, ap; options...)
    L = samplingparameter(sstyle, ap; options...)
    oversampling_grid(ap, L; options...)
end

sampling_grid(samplingstyle::ProductSamplingStyle, ap; options...) =
    tensorproduct( map((x,style) -> sampling_grid(x; samplingstyle=style, options...),
                        elements(ap), samplingstyle.styles)...
    )


## Discretization

# The discretization requires the sampling operator. If it is not supplied, we
# compute it first.
discretization(ap::ApproximationProblem; options...) =
    discretization(ap, samplingoperator(ap; options...); options...)

discretization(ap::ApproximationProblem, S; options...) = apply(S, dictionary(ap); options...)

# Optionally, supplying f as the first argument yields the right hand side as well
discretization(f, ap::ApproximationProblem; options...) =
    discretization(f, ap, samplingoperator(ap; options...); options...)

function discretization(f, ap::ApproximationProblem, S; options...)
    A = discretization(ap, S; options...)
    B = apply(S, f; options...)
    A, B
end


## Dual dictionary and sampling operator

dualdictionary(ap::DictionaryApproximation; options...) = dualdictionary(dictionary(ap))
dualdictionary(ap::PlatformApproximation; options...) =
    dualdictionary(ap.platform, ap.param; dict = dictionary(ap))

discrete_normalization(ap::DictionaryApproximation; options...) =
    discrete_normalization(dictionary(ap), samplingparameter(ap); options...)
discrete_normalization(ap::PlatformApproximation; options...) =
    discrete_normalization(ap.platform, ap.param, samplingparameter(ap); options...)


dualdiscretization(ap::ApproximationProblem; options...) =
    dualdiscretization(ap, dualsamplingoperator(ap; options...); options...)

dualdiscretization(ap::ApproximationProblem, Stilde; options...) =
    apply(Stilde, dualdictionary(ap; options...); options...)


##  Solver

function solver(ap::ApproximationProblem;
            samplingstyle = SamplingStyle(ap),
            solverstyle = SolverStyle(ap, samplingstyle), options...)
    S = samplingoperator(samplingstyle, ap; options...)
    A = discretization(ap, S; options...)
    solver(solverstyle, ap, A; S = S, options...)
end



## The AZ algorithm

AZ_A(ap::ApproximationProblem; options...) = discretization(ap; options...)

AZ_Z(ap::ApproximationProblem; options...) = dualdiscretization(ap; options...)
AZ_Zt(ap::ApproximationProblem; options...) = AZ_Z(ap; options...)'

function plungeoperator(ap::ApproximationProblem; options...)
    A = discretization(ap; options...)
    Z = dualdiscretization(ap; options...)
    I = IdentityOperator(dest(A))
    I - A*Z'
end

function plungematrix(ap::ApproximationProblem; options...)
    A = discretization(ap; options...)
    P = plungeoperator(ap; options...)
    P * A
end

function plungerank(ap::ApproximationProblem;
            REG = default_regularization,
            rankestimate = 40,
            threshold = regularization_threshold(coefficienttype(ap)), options...)
    C = plungematrix(ap; options...)
    Q = REG(C; threshold = threshold, rankestimate = rankestimate, options...)
    length(Q.Sinv)
end


# # TODO: clean up these scaling factors. They have to do with a normalization of the
# # sampling operators.
Zt_scaling_factor(S::Dictionary, A) = length(supergrid(grid(dest(A))))
Zt_scaling_factor(S::DerivedDict, A) = Zt_scaling_factor(superdict(S), A)
Zt_scaling_factor(S::ChebyshevT, A) = length(supergrid(grid(dest(A))))/2

smoothingoperator(ap::ApproximationProblem; options...) =
    WeightedSmoothingOperator(dictionary(ap); options...)
