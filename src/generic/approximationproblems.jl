##########################
# Approximation problems
##########################

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


# Caching of properties
cache(ap::ApproximationProblem) = ap.cache
cache(ap::ApproximationProblem, property) = ap.cache[property]
cache!(ap::ApproximationProblem, property, value) = ap.cache[property] = value

ap_hasproperty(ap::ApproximationProblem, property) = haskey(cache(ap), property)

ap_fun(ap::ApproximationProblem) = cache(ap, ap_fun)
ap_fun!(ap::ApproximationProblem, f) = cache!(ap, ap_fun, f)

# For convenience
measure(ap::ApproximationProblem) = measure(platform(ap))
coefficienttype(ap::ApproximationProblem) = coefficienttype(platform(ap))

approximation_measure(ap::ApproximationProblem) = approximation_measure(SamplingStyle(ap), ap)
approximation_measure(ss::DiscreteStyle, ap::ApproximationProblem) = discretemeasure(ap)
approximation_measure(ss::ProjectionStyle, ap::ApproximationProblem) = measure(ap)
approximation_measure(ss::ProductSamplingStyle, ap::ApproximationProblem) = approximation_measure(factor(ss,1),ap)



"A `PlatformApproximation` stores a platform and a platform parameter."
mutable struct PlatformApproximation <: ApproximationProblem
    platform    ::  Platform
    plt_par         # platform parameter
    cache

    PlatformApproximation(platform, plt_par) = new(platform, plt_par, Dict{Any,Any}())
    function PlatformApproximation(platform, plt_par, smpl_par)
        cache = Dict{Any,Any}()
        cache[samplingparameter] = smpl_par
        new(platform, plt_par, cache)
    end
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

function dualdictionary(ap::PlatformApproximation)
    if ap_hasproperty(ap, dualdictionary)
        cache(ap, dualdictionary)
    else
        Φ = dualdictionary(platform(ap), platformparameter(ap), approximation_measure(ap))
        cache!(ap, dualdictionary, Φ)
    end
end

coefficienttype(ap::PlatformApproximation) = coefficienttype(dictionary(ap))


# specialized functions for when the platform is a product platform

productparameter(ap::PlatformApproximation) = productparameter(platform(ap), platformparameter(ap))

function factors(ap::PlatformApproximation)
    if ap_hasproperty(ap, samplingparameter)
        smpl_par = samplingparameter(ap)
        if length(smpl_par)==length(productparameter(ap))
            map(approximationproblem, factors(platform(ap)), productparameter(ap), smpl_par)
        else
            error("sampling parameter should contain $(length(productparameter(ap))) elements")
        end
    else
        map(approximationproblem, components(platform(ap)), productparameter(ap))
    end
end

unsafe_factors(ap::PlatformApproximation) =
    map(approximationproblem, components(platform(ap)), productparameter(ap))

components(ap::PlatformApproximation) = factors(ap)


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
to_platform(dict::Dictionary, smpl_par) =
    (platform(dict), platformparameter(dict), smpl_par)
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

function copy_cache(src, dest)
    for k in keys(cache(src))
        cache!(dest, k, cache(src, k))
    end
end

# From the adaptive approximation problem we can create a concrete platform
# approximation by specifying a parameter value.
function approximationproblem(ap::AdaptiveApproximation, plt_par)
    ap2 = approximationproblem(platform(ap), plt_par)
    copy_cache(ap, ap2)
    ap2
end
function approximationproblem(ap::AdaptiveApproximation, plt_par, smpl_par)
    ap2 = approximationproblem(platform(ap), plt_par, smpl_par)
    copy_cache(ap, ap2)
    ap2
end


# 4. Enable options
function approximationproblem(args...; options...)
    ap = approximationproblem(args...)
    ap_process_options(ap; options...)
    ap
end

# 5. Allow a function to be specified as first argument
approximationproblem(fun, dict::Dictionary, args...) =
    approximationproblem(fun, approximationproblem(dict, args...))
approximationproblem(fun, platform::Platform, args...) =
    approximationproblem(fun, approximationproblem(platform, args...))
function approximationproblem(fun, ap::ApproximationProblem)
    ap_fun!(ap, fun)
    ap
end

export _ap
const _ap = approximationproblem
