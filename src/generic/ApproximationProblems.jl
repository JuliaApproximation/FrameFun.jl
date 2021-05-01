##########################
# Approximation problems
##########################
module ApproximationProblems
using ..Platforms
import ..Platforms: SamplingStyle, SolverStyle, dictionary, platform
import BasisFunctions: coefficienttype, component, components
using BasisFunctions: Dictionary, TensorProductDict
using ..ExtensionFrames
using DomainSets: Domain

export ApproximationProblem
"""
    abstract type ApproximationProblem end

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

export samplingparam
"""
    samplingparam(ap::ApproximationProblem)

Return the sampling parameter if it was saved by the approximation problem.
Otherwise return `nothing`
"""
samplingparam(ap::ApproximationProblem) = ap.samplingparam

export setsamplingparam!
"""
    setsamplingparam!(ap::ApproximationProblem, L)

Save the sampling parameter in the approximation problem
"""
setsamplingparam!(ap::ApproximationProblem, L) = (ap.samplingparam = L)

export AbstractPlatformApproximation

"""
    abstract type AbstractPlatformApproximation <: ApproximationProblem

An approximation problem with a `platform` and a `parameter`
"""
abstract type AbstractPlatformApproximation <: ApproximationProblem end
platform(ap::AbstractPlatformApproximation) = ap.platform
export parameter
parameter(ap::AbstractPlatformApproximation) = ap.param
coefficienttype(ap::AbstractPlatformApproximation) = coefficienttype(ap.dict)

export PlatformApproximation
"""
    mutable struct PlatformApproximation <: ApproximationProblem

A `PlatformApproximation` stores a platform and a parameter value. This
corresponds to a specific dictionary (which is also stored). All queries are
delegated to the platform, and thus the answers may differ from the defaults.

As with dictionary approximations, the settings can be overruled by explicitly
specifying optional arguments to `Fun`.
"""
mutable struct PlatformApproximation <: AbstractPlatformApproximation
    platform    ::  Platform
    dict        ::  Dictionary
    param
    samplingparam

    PlatformApproximation(platform, param, samplingparam = nothing) =
        new(platform, dictionary(platform, param), param, samplingparam)
end

dictionary(ap::PlatformApproximation) = ap.dict


SamplingStyle(ap::PlatformApproximation) = SamplingStyle(ap.platform)
SolverStyle(samplingstyle::SamplingStyle, ap::PlatformApproximation) = SolverStyle(ap.platform, samplingstyle)


export AdaptiveApproximation
"""
    struct AdaptiveApproximation <: ApproximationProblem

An `AdaptiveApproximation` only stores a platform, with which to compute adaptive
approximations.
"""
struct AdaptiveApproximation <: ApproximationProblem
    platform    ::  Platform
end

dictionary(ap::AdaptiveApproximation, param) = ap.platform[param]

export DictionaryApproximation
mutable struct DictionaryApproximation <: ApproximationProblem
end


export ProductPlatformApproximation
"""
    mutable struct ProductPlatformApproximation{N} <: ApproximationProblem

A `ProductPlatformApproximation` corresponds to a product platform.
"""
mutable struct ProductPlatformApproximation{N} <: AbstractPlatformApproximation
    platform        ::  ProductPlatform{N}
    param
    productparam    ::  NTuple{N,Any}
    dict            ::  TensorProductDict

    samplingparam

    ProductPlatformApproximation{N}(platform::ProductPlatform, param, samplingparam=nothing) where {N} =
        new(platform, param, productparameter(platform, param), dictionary(platform, param), samplingparam)
end

setsamplingparam!(ap::ProductPlatformApproximation, param) =
    (@assert length(ap.productparam)==length(param);ap.samplingparam = param)

approximationproblem(platform::ProductPlatform{N}, param) where {N} =
    ProductPlatformApproximation{N}(platform, param)

SamplingStyle(ap::ProductPlatformApproximation) = SamplingStyle(ap.platform)
SolverStyle(samplingstyle::SamplingStyle, ap::ProductPlatformApproximation) = SolverStyle(ap.platform, samplingstyle)

components(ap::ProductPlatformApproximation) = nothing == samplingparam(ap) ?
    error("Not possible to get components of `ProductPlatformApproximation` without `samplingparam`") :
    length(samplingparam(ap))==length(ap.productparam) ?
    map(approximationproblem, components(ap.platform), ap.productparam, samplingparam(ap)) :
    error("sampling parameter should contain $(length(ap.productparam)) elements")


unsafe_components(ap::ProductPlatformApproximation) = map(approximationproblem, components(ap.platform), ap.productparam)

dictionary(ap::ProductPlatformApproximation) = ap.dict


export approximationproblem
"""
    approximationproblem(args...)

# examples

## Create an approxmationproblem from a dictionary
```jldocs
julia> approximationproblem(Fourier(10))
```

## Create an approxmationproblem from a dictionary and a domain (extensionframe)
```jldocs
julia> approximationproblem(Fourier(10),0.0..0.5)
```

## Create an approxmationproblem from a platform
```jldocs
julia> approximationproblem(FourierPlatform())
```

## Create an approxmationproblem from a platform and a platform parameter
```jldocs
julia> approximationproblem(FourierPlatform(), 10)
```
"""
approximationproblem(dict::Dictionary) =
    approximationproblem(platform(dict), param(dict))
# This two-argument routine allows to specify a coefficient type, in which case the
# dictionary will be promoted if necessary. This allows the `Fun` constructor to ensure
# that the dictionary can handle complex coefficients, for example.
# approximationproblem(::Type{T}, dict::Dictionary) where {T} =
#     DictionaryApproximation(ensure_coefficienttype(T, dict))



# 2. An approximation problem from a platform and a concrete parameter value
approximationproblem(platform::Platform, param) = (@assert correctparamformat(platform, param);PlatformApproximation(platform, param))
approximationproblem(platform::Platform, param, L) = (@assert correctparamformat(platform, param);PlatformApproximation(platform, param, L))

# 3. An adaptive approximation problem from a platform
approximationproblem(platform::Platform) = AdaptiveApproximation(platform)

# From the adaptive approximation problem we can create a concrete platform
# approximation by specifying a parameter value.
approximationproblem(ap::AdaptiveApproximation, param) = approximationproblem(ap.platform, param)
approximationproblem(ap::AdaptiveApproximation, param, L) = approximationproblem(ap.platform, param, L)
end
