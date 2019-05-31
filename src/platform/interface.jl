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

elements(ap::PlatformApproximation) = map(PlatformApproximation, elements(platform(ap)), elements(parameter(ap)))

element(ap::PlatformApproximation, i) = PlatformApproximation(element(platform(ap), i), element(parameter(ap), i))

elements(t::NTuple) = t

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
approximationproblem(platform::Platform, param) = PlatformApproximation(platform, param)
approximationproblem(platform::Platform, param, L) = PlatformApproximation(platform, param, L)

# 3. An adaptive approximation problem from a platform
approximationproblem(platform::Platform) = AdaptiveApproximation(platform)

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
for op in (:discretization, :dualdiscretization, :measure, :dualplatformdictionary)
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

function interpolation_grid(platform::Platform, n; dict = dictionary(platform, n), oversamplingfactor=1, options...)
    oversamplingfactor != 1 && @warn "InterpolationStyle does not support option `oversamplingfactor=$oversamplingfactor`, use `OversamplingStyle()` instead."
    interpolation_grid(dict)
end

oversampling_grid(platform::Platform, n, L; dict = dictionary(platform, n), options...) =
    oversampling_grid(dict, L)

space(platform::Platform) = space(measure(platform))


########################
# Interface delegation
########################

# Here, we direct queries to the dictionary or to the platform.
# Also, queries to the platform are then redirected to the dictionary. Some
# of the redirection is implement in platform.jl and not here.

## SamplingStyle and SolverStyle

SamplingStyle(ap::DictionaryApproximation) = SamplingStyle(dictionary(ap))
SolverStyle(samplingstyle::SamplingStyle, ap::DictionaryApproximation) = SolverStyle(dictionary(ap), samplingstyle)

SamplingStyle(ap::PlatformApproximation) = SamplingStyle(ap.platform)
SolverStyle(samplingstyle::SamplingStyle, ap::PlatformApproximation) = SolverStyle(ap.platform, samplingstyle)

## The sampling parameter

samplingparameter(ap::ApproximationProblem; samplingstyle = SamplingStyle(ap), options...) =
    samplingparameter(samplingstyle, ap; options...)

samplingparameter(samplingstyle::DiscreteGramStyle, ap::ApproximationProblem; options...) =
    samplingparameter(OversamplingStyle(), ap; options...)

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
            verbose = false, oversamplingfactor = 1, options...)
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
            T = coefficienttype(ap), options...)
            # grid = sampling_grid(samplingstyle, ap; options...)
    dmeasure = haskey(options,:samplingmeasure) ? options[:samplingmeasure] :
                    discretemeasure(samplingstyle, ap; options...)
    g  = grid(dmeasure)
    w = BasisFunctions.weights(dmeasure)
    GS = GridSampling(g, T)
    D = DiagonalOperator(dest(GS), dest(GS), sqrt.(w))
    D*GS
end

samplingoperator(samplingstyle::DiscreteGramStyle, ap::ApproximationProblem; options...) =
    ProjectionSampling(dictionary(ap), measure(samplingstyle, ap; options...))

samplingoperator(::GenericSamplingStyle, ap::PlatformApproximation; options...) =
    genericsamplingoperator(ap.platform, ap.param; dict = ap.dict, options...)

samplingoperator(samplingstyle::GramStyle, ap::ApproximationProblem; options...) =
    ProjectionSampling(dictionary(ap), measure(samplingstyle, ap; options...))

dualsamplingoperator(samplingstyle::GramStyle, ap::ApproximationProblem; options...) =
    ProjectionSampling(
        haskey(options, :dualdict) ? options[:dualdict] : dualplatformdictionary(samplingstyle, ap; options...),
            measure(samplingstyle, ap; options...))

dualsamplingoperator(samplingstyle::DiscreteGramStyle, ap::ApproximationProblem; options...) =
    ProjectionSampling(
        haskey(options, :dualdict) ? options[:dualdict] : dualplatformdictionary(samplingstyle, ap; options...),
            measure(samplingstyle, ap; options...))

samplingoperator(samplingstyle::RectangularGramStyle, ap::ApproximationProblem;
            projectiondict, options...) =
    ProjectionSampling(projectiondict, measure(samplingstyle, ap; options...))

samplingoperator(samplingstyle::ProductSamplingStyle, ap::ApproximationProblem; options...) =
    tensorproduct( map( (x,style) -> samplingoperator(x; samplingstyle=style, options...),
                        elements(ap), samplingstyle.styles)...
    )

# First, if necessary, we recompute the samplingoperator using the given options.
function dualsamplingoperator(samplingstyle::DiscreteStyle, ap::ApproximationProblem; options...)
    if haskey(options, :S)
        S = options[:S]
    else
        S = samplingoperator(samplingstyle, ap; options...)
    end
    S
end

dualsamplingoperator(samplingstyle::SamplingStyle, ap::ApproximationProblem; options...) =
    samplingoperator(samplingstyle, ap; options...)

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
function sampling_grid(::InterpolationStyle, ap::DictionaryApproximation; oversamplingfactor=1, options...)
    oversamplingfactor != 1 && @warn "InterpolationStyle does not support option `oversamplingfactor=$oversamplingfactor`, use `OversamplingStyle()` instead."
    interpolation_grid(dictionary(ap))
end
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

function sampling_grid(sstyle::OversamplingStyle, ap::ApproximationProblem; options...)
    L = samplingparameter(sstyle, ap; options...)
    oversampling_grid(ap, L; options...)
end
sampling_grid(sstyle::DiscreteGramStyle, ap; options...) =
    sampling_grid(OversamplingStyle(), ap; options...)

sampling_grid(samplingstyle::ProductSamplingStyle, ap; options...) =
    tensorproduct( map((x,style) -> sampling_grid(x; samplingstyle=style, options...),
                        elements(ap), samplingstyle.styles)...
    )


## Discretization
# The discretization requires the sampling operator only in the discrete case.
discretization(ap::ApproximationProblem; samplingstyle = SamplingStyle(ap), options...) =
    discretization(samplingstyle, ap; options...)

# The discretization requires the sampling operator. If it is not supplied, we
# compute it first.
discretization(sstyle::SamplingStyle, ap::ApproximationProblem; options...) =
    discretization(sstyle, ap, samplingoperator(sstyle, ap; options...); options...)

discretization(::SamplingStyle, ap::ApproximationProblem, S; options...) =
    apply(S, dictionary(ap); options...)


# Optionally, supplying f as the first argument yields the right hand side as well
discretization(f, ap::ApproximationProblem; samplingstyle=SamplingStyle(ap), options...) =
    discretization(f, samplingstyle, ap; options...)

discretization(f, sstyle::SamplingStyle, ap; options...) =
    discretization(f, sstyle, ap, samplingoperator(sstyle, ap; options...); options...)

function discretization(f, sstyle::SamplingStyle, ap::ApproximationProblem, S; options...)
    A = discretization(sstyle, ap, S; options...)
    B = apply(S, f; options...)
    A, B
end

dualdiscretization(ap::ApproximationProblem; samplingstyle=SamplingStyle(ap), options...) =
    dualdiscretization(samplingstyle, ap; options...)

function dualdiscretization(sstyle::SamplingStyle, ap::ApproximationProblem; options...)
    dualdict = haskey(options, :dualdict) ? options[:dualdict] : dualplatformdictionary(sstyle, ap; options...)
    dualdiscretization(sstyle, ap, dualsamplingoperator(sstyle, ap; dualdict=dualdict, options...); dualdict=dualdict, options...)
end

function dualdiscretization(sstyle::SamplingStyle, ap::ApproximationProblem, Stilde; options...)
    dualdict = haskey(options, :dualdict) ? options[:dualdict] : dualplatformdictionary(sstyle, ap; options...)
    apply(Stilde, dualdict; options...)
end

normalizationoperator(::DictionaryOperatorStyle, ap::ApproximationProblem; sstyle=SamplingStyle(ap), opts...) =
    normalizationoperator(sstyle, ap; opts...)

function normalizationoperator(sstyle::DiscreteStyle, ap::ApproximationProblem; T=coefficienttype(ap), options...)
    grid = sampling_grid(sstyle, ap; options...)
    sampling_normalization(GridBasis{T}(grid),
        measure(sstyle, ap; options...))
end

##  Solver

solver(ap::ApproximationProblem;
            problemstyle = ProblemStyle(ap),
            samplingstyle = SamplingStyle(ap),
            solverstyle = SolverStyle(samplingstyle, ap), options...) =
        solver(problemstyle, solverstyle, ap; samplingstyle=samplingstyle, options...)

function solver(pstyle::DictionaryOperatorStyle, solverstyle::SolverStyle, ap::ApproximationProblem;
            samplingstyle = SamplingStyle(ap), options...)
    S = samplingoperator(samplingstyle, ap; options...)
    A = discretization(samplingstyle, ap, S; options...)
    solver(solverstyle, ap, A; S = S, options...)
end


## The AZ algorithm
AZ_A(ap::ApproximationProblem; problemstyle=ProblemStyle(ap), options...) =
    AZ_A(problemstyle, ap; options...)
AZ_Z(ap::ApproximationProblem; problemstyle=ProblemStyle(ap), options...) =
    AZ_Z(problemstyle, ap; options...)
AZ_Zt(ap::ApproximationProblem; problemstyle=ProblemStyle(ap), options...) =
    AZ_Zt(problemstyle, ap; options...)
plungeoperator(ap::ApproximationProblem; problemstyle=ProblemStyle(ap), options...) =
    plungeoperator(problemstyle, ap; options...)
plungematrix(ap::ApproximationProblem; problemstyle=ProblemStyle(ap), options...) =
    plungematrix(problemstyle, ap; options...)

default_aznormalization(a...) = false

AZ_A(pstyle::ProblemStyle, ap; normalizedsampling=default_aznormalization(ap), options...) =
    normalizedsampling ?
        normalizationoperator(pstyle, ap;options...)*_AZ_A(pstyle, ap; options...) :
        _AZ_A(pstyle, ap; options...)

AZ_Z(pstyle::ProblemStyle, ap; normalizedsampling=default_aznormalization(ap), options...) =
    normalizedsampling ?
        inv(normalizationoperator(pstyle,ap; options...))*_AZ_Z(pstyle, ap; options...) :
        _AZ_Z(pstyle, ap; options...)

AZ_Zt(pstyle::ProblemStyle, ap::ApproximationProblem; options...) = AZ_Z(pstyle, ap; options...)'

_AZ_A(::DictionaryOperatorStyle, ap; samplingstyle = SamplingStyle(ap), options...) =
    discretization(samplingstyle, ap, samplingoperator(samplingstyle, ap; options...); options...)
_AZ_Z(::DictionaryOperatorStyle, ap; samplingstyle = SamplingStyle(ap), options...) =
    dualdiscretization(samplingstyle, ap, dualsamplingoperator(samplingstyle, ap; options...); options...)

function plungeoperator(problemstyle::ProblemStyle, ap::ApproximationProblem; options...)
    A = AZ_A(problemstyle, ap; options...)
    Zt = AZ_Zt(problemstyle, ap; options...)
    # A = discretization(ap; options...)
    # Z = dualdiscretization(ap; options...)
    # I = IdentityOperator(dest(A))
    I - A*Zt
end


function plungematrix(pstyle::ProblemStyle, ap::ApproximationProblem; options...)
    A = AZ_A(pstyle, ap; options...)
    P = plungeoperator(pstyle, ap; options...)
    P * A
end

# function plungematrix(::DictionaryOperatorStyle, ap::ApproximationProblem; options...)
#     A = discretization(ap; options...)
#     P = plungeoperator(ap; options...)
#     P * A
# end

function plungerank(ap::ApproximationProblem;
            REG = default_regularization,
            rankestimate = 40,
            threshold = regularization_threshold(coefficienttype(ap)), options...)
    C = plungematrix(ap; options...)
    Q = REG(C; threshold = threshold, rankestimate = rankestimate, options...)
    length(Q.Sinv)
end

smoothingoperator(ap::ApproximationProblem; options...) =
    WeightedSmoothingOperator(dictionary(ap); options...)

###
# GenericOperator Interface
###

## The AZ algorithm
AZ_A(::GenericOperatorStyle, ap::ApproximationProblem; options...) = SynthesisOperator(dictionary(ap), measure(ap; options...))
AZ_Z(::GenericOperatorStyle, ap::ApproximationProblem; options...) = SynthesisOperator(dualplatformdictionary(ap; options...), measure(ap; options...))
samplingoperator(pstyle::GenericOperatorStyle, ap; options...) = AZ_Zt(pstyle, ap; options...)
samplingoperator(pstyle::DictionaryOperatorStyle, ap; options...) = samplingoperator(ap; options...)

function plungematrix(pstyle::GenericOperatorStyle, ap::ApproximationProblem; options...)
    S = samplingoperator(pstyle, ap; options...)
    A = AZ_A(pstyle, ap; options...)
    Zt = AZ_Zt(pstyle, ap; options...)
    tmp1 = apply(Zt, A; options...)
    tmp2 = apply(S, A; options...)
    tmp2 - tmp2*tmp1
end
