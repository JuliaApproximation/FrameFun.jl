##########################
# Approximation problems
##########################

import BasisFunctions: coefficienttype
using BasisFunctions: Dictionary, TensorProductDict
using DomainSets: Domain

"""
    abstract type ApproximationProblem end

Approximation problem types group the arguments that a user supplies
to the `Fun` constructor. This may be a dictionary or a platform, for example.
By passing around approximation problem objects, functions can easily
implement the same interface as `Fun`.

All functionality is eventually delegated back to the underlying dictionary or
platform. Thus, an approximation problem is nothing but an empty intermediate
layer that passes on the questions to the right object regardless of how the
question was asked.
"""
abstract type ApproximationProblem end

"Return the sampling parameter if it was set, otherwise return `nothing`."
_samplingparameter(ap::ApproximationProblem) = ap.smpl_par

"Set the sampling parameter in the approximation problem."
_samplingparameter!(ap::ApproximationProblem, smpl_par) = (ap.smpl_par = smpl_par)


"Supertype of approximation problems based on a `platform`."
abstract type PlatformApproximationProblem <: ApproximationProblem end

platform(ap::PlatformApproximationProblem) = ap.platform

export platformparameter
"Return the platform parameter of the approximation problem."
platformparameter(ap::PlatformApproximationProblem) = ap.plt_par

coefficienttype(ap::PlatformApproximationProblem) = coefficienttype(dictionary(ap))

"""
    mutable struct PlatformApproximation <: PlatformApproximationProblem

A `PlatformApproximation` stores a platform and a platform parameter. This
corresponds to a specific dictionary (which is also stored). All queries are
delegated to the platform, and thus the answers may differ from the defaults.

As with dictionary approximations, the settings can be overruled by explicitly
specifying optional arguments to `Fun`.
"""
mutable struct PlatformApproximation <: PlatformApproximationProblem
    platform    ::  Platform
    dict        ::  Dictionary
    plt_par
    smpl_par
    cache

    PlatformApproximation(platform, plt_par, smpl_par = nothing) =
        new(platform, dictionary(platform, plt_par), plt_par, smpl_par)
end

dictionary(ap::PlatformApproximation) = ap.dict


SamplingStyle(ap::PlatformApproximation) = SamplingStyle(ap.platform)
SolverStyle(samplingstyle::SamplingStyle, ap::PlatformApproximation) =
    SolverStyle(ap.platform, samplingstyle)


"""
    struct AdaptiveApproximation <: ApproximationProblem

An `AdaptiveApproximation` only stores a platform, with which to compute adaptive
approximations. It has no single platform parameter.
"""
struct AdaptiveApproximation <: ApproximationProblem
    platform    ::  Platform
    cache
end

platform(ap::AdaptiveApproximation) = ap.platform

dictionary(ap::AdaptiveApproximation, plt_par) = dictionary(platform(ap), plt_par)

platformparameter(ap::AdaptiveApproximation) =
    error("An adaptive approximation problem does not contain a platform parameter.")


"""
    mutable struct ProductPlatformApproximation{N} <: ApproximationProblem

A `ProductPlatformApproximation` corresponds to a product platform.
"""
mutable struct ProductPlatformApproximation{N} <: PlatformApproximationProblem
    platform        ::  ProductPlatform{N}
    plt_par
    productparam    ::  NTuple{N,Any}
    dict            ::  TensorProductDict
    smpl_par
    cache


    function ProductPlatformApproximation{N}(platform::ProductPlatform, plt_par, smpl_par=nothing) where {N}
        new(platform, plt_par, productparameter(platform, plt_par), dictionary(platform, plt_par), smpl_par)
    end
end

function _samplingparameter!(ap::ProductPlatformApproximation, smpl_par)
    @assert length(ap.productparam)==length(smpl_par)
    ap.smpl_par = smpl_par
end

approximationproblem(platform::ProductPlatform{N}, plt_par) where {N} =
    ProductPlatformApproximation{N}(platform, plt_par)

SamplingStyle(ap::ProductPlatformApproximation) = SamplingStyle(ap.platform)
SolverStyle(samplingstyle::SamplingStyle, ap::ProductPlatformApproximation) =
    SolverStyle(ap.platform, samplingstyle)

function components(ap::ProductPlatformApproximation)
    smpl_par = _samplingparameter(ap)
    if smpl_par == nothing
        error("Not possible to get components of `ProductPlatformApproximation` without `samplingparam`")
    end
    if length(smpl_par)==length(ap.productparam)
        map(approximationproblem, components(ap.platform), ap.productparam, smpl_par)
    else
        error("sampling parameter should contain $(length(ap.productparam)) elements")
    end
end


unsafe_components(ap::ProductPlatformApproximation) =
    map(approximationproblem, components(ap.platform), ap.productparam)

dictionary(ap::ProductPlatformApproximation) = ap.dict



to_platform(dict::Dictionary) = (platform(dict), platformparameter(dict))
function to_platform(dict::Dictionary, plt_par)
    plt = platform(dict)
    @assert correctparamformat(plt, plt_par)
    plt, plt_par
end
function to_platform(dict::Dictionary, plt_par, smpl_par)
    plt = platform(dict)
    @assert correctparamformat(plt, plt_par)
    plt, plt_par, smpl_par
end
function to_platform(dict::Dictionary, domain::Domain, args...)
    if domain == support(dict)
        to_platform(dict, args...)
    else
        to_platform(extensionframe(domain, dict), args...)
    end
end


export approximationproblem
"""
    approximationproblem(args...)

# examples

## Create an approximationproblem from a dictionary
```jldocs
julia> approximationproblem(Fourier(10))
```

## Create an approximationproblem from a dictionary and a domain (extensionframe)
```jldocs
julia> approximationproblem(Fourier(10),0.0..0.5)
```

## Create an approximationproblem from a platform
```jldocs
julia> approximationproblem(FourierPlatform())
```

## Create an approximationproblem from a platform and a platform parameter
```jldocs
julia> approximationproblem(FourierPlatform(), 10)
```
"""
approximationproblem(dict::Dictionary, args...) =
    approximationproblem(to_platform(dict, args...)...)


# 2. An approximation problem from a platform and a concrete parameter value
function approximationproblem(platform::Platform, plt_par)
    @assert correctparamformat(platform, plt_par)
    PlatformApproximation(platform, plt_par)
end
function approximationproblem(platform::Platform, plt_par, smpl_par)
    @assert correctparamformat(platform, plt_par)
    PlatformApproximation(platform, plt_par, smpl_par)
end

# 3. An adaptive approximation problem from a platform
approximationproblem(platform::Platform) = AdaptiveApproximation(platform)

# From the adaptive approximation problem we can create a concrete platform
# approximation by specifying a parameter value.
approximationproblem(ap::AdaptiveApproximation, plt_par) =
    approximationproblem(ap.platform, plt_par)
approximationproblem(ap::AdaptiveApproximation, plt_par, smpl_par) =
    approximationproblem(ap.platform, plt_par, smpl_par)
