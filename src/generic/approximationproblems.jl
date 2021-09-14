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

Some of the properties of an approximation problem may be cached in order to
improve efficiency.
"""
abstract type ApproximationProblem end


"Return the sampling parameter if it was set, otherwise return `nothing`."
_samplingparameter(ap::ApproximationProblem) = ap.smpl_par

"Set the sampling parameter in the approximation problem."
function _samplingparameter!(ap::ApproximationProblem, smpl_par)
    if platform(ap) isa ProductPlatform
        @assert length(productparameter(ap))==length(smpl_par)
    end
    (ap.smpl_par = smpl_par)
end


SamplingStyle(ap::ApproximationProblem) = SamplingStyle(platform(ap))
SolverStyle(samplingstyle::SamplingStyle, ap::ApproximationProblem) =
    SolverStyle(platform(ap), samplingstyle)


# Caching of properties
cache(ap::ApproximationProblem) = ap.cache
cache(ap::ApproximationProblem, property) = ap.cache[property]
cache!(ap::ApproximationProblem, property, value) = ap.cache[property] = value

ap_hasproperty(ap::ApproximationProblem, property) = haskey(cache(ap), property)

ap_fun(ap::ApproximationProblem) = cache(ap, ap_fun)
ap_fun!(ap::ApproximationProblem, f) = cache!(ap, ap_fun, f)



"A `PlatformApproximation` stores a platform and a platform parameter."
mutable struct PlatformApproximation <: ApproximationProblem
    platform    ::  Platform
    plt_par         # platform parameter
    smpl_par        # sampling parameter
    cache

    PlatformApproximation(platform, plt_par, smpl_par = nothing) =
        new(platform, plt_par, smpl_par, Dict{Any,Any}())
end

platform(ap::PlatformApproximation) = ap.platform

export platformparameter
"Return the platform parameter of the approximation problem."
platformparameter(ap::PlatformApproximation) = ap.plt_par

function dictionary(ap::PlatformApproximation)
    if ap_hasproperty(ap, dictionary)
        cache(ap, dictionary)
    else
        Φ = dictionary(platform(ap), platformparameter(ap))
        cache!(ap, dictionary, Φ)
    end
end


coefficienttype(ap::PlatformApproximation) = coefficienttype(dictionary(ap))

# specialized functions for when the platform is a product platform

productparameter(ap::PlatformApproximation) = productparameter(platform(ap), platformparameter(ap))

function productcomponents(ap::PlatformApproximation)
    smpl_par = _samplingparameter(ap)
    if smpl_par == nothing
        error("Not possible to get product components of approximation problem without known sampling parameter")
    end
    if length(smpl_par)==length(productparameter(ap))
        map(approximationproblem, components(ap.platform), productparameter(ap), smpl_par)
    else
        error("sampling parameter should contain $(length(productparameter(ap))) elements")
    end
end

unsafe_productcomponents(ap::PlatformApproximation) =
    map(approximationproblem, components(platform(ap)), productparameter(ap))



"""
    struct AdaptiveApproximation <: ApproximationProblem

An `AdaptiveApproximation` only stores a platform, with which to compute adaptive
approximations. It has no single platform parameter.
"""
struct AdaptiveApproximation <: ApproximationProblem
    platform    ::  Platform
    cache
    AdaptiveApproximation(platform::Platform) = new(platform, Dict{Any,Any}())
end

platform(ap::AdaptiveApproximation) = ap.platform

dictionary(ap::AdaptiveApproximation) =
    error("An adaptive approximation problem does not have a single dictionary.")
dictionary(ap::AdaptiveApproximation, plt_par) = dictionary(platform(ap), plt_par)

platformparameter(ap::AdaptiveApproximation) =
    error("An adaptive approximation problem does not contain a platform parameter.")



## Now on to the actual interface:

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
    approximationproblem(platform(ap), plt_par)
approximationproblem(ap::AdaptiveApproximation, plt_par, smpl_par) =
    approximationproblem(platform(ap), plt_par, smpl_par)


# Allow a function to be specified as first argument
approximationproblem(fun, dict::Dictionary, args...) =
    approximationproblem(fun, to_platform(dict, args...)...)
function approximationproblem(fun, platform::Platform, args...)
    ap = approximationproblem(platform, args...)
    ap_fun!(ap, fun)
    ap
end
