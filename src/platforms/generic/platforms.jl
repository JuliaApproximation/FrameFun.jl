
using BasisFunctions: Dictionary, TensorProductDict, Dictionary1d, tensorproduct,
    extensionsize, gramdual, Weight, ProductWeight, DiscreteProductWeight, resize,
    isbasis, hastransform, dimensions, productmeasure
import Base: getindex
import BasisFunctions: components, factors, dictionary, measure, Measure

#################
# Platform
#################

export Platform, BasisPlatform, FramePlatform
"""
    abstract type Platform end

A platform represents a family of dictionaries.

A platform typically has a primal and a dual sequence. The platform maps an
index to a set of parameter values, which is then used to generate a primal
or dual dictionary (depending on the chosen measure).

See also [`BasisPlatform`](@ref), `FramePlatform`](@ref)
"""
abstract type Platform end

getindex(platform::Platform, param) = dictionary(platform, param)

correctparamformat(platform::Platform, _) = false

function dictionary(platform::Platform, param)
    if !correctparamformat(platform, param)
        throw(ArgumentError("Parameter $param not suited for platform $platform. "))
    end
    unsafe_dictionary(platform, param)
end

export dictionary, dualdictionary

"""
    dualdictionary(platform::Platform, param, measure::Measure; options...)

Return the dual dictionary of the platform.
"""
dualdictionary(platform::Platform, param, measure::Measure; dict=dictionary(platform, param), options...) =
    gramdual(dict, measure; options...)

include("platformadaptivity.jl")

"""
    abstract type BasisPlatform <: Platform end

A `BasisPlatform` represents a family of bases.
"""
abstract type BasisPlatform <: Platform end

"""
    abstract type FramePlatform <: Platform end

A `FramePlatform` represents a family of frames
(or ill-conditioned bases which are in the limit a frame).
"""
abstract type FramePlatform <: Platform end


include("styles.jl")




#################
# ModelPlatform
#################

"""
    struct ModelPlatform <: Platform

A `ModelPlatform` is a platform based on a model dictionary. The platform is
defined by resizing the dictionary, using its own implementation of `resize`.
All other operations are the defaults for the model dictionary.

This platform is convenient to compute adaptive approximations based on an
example of a dictionary from the desired family.
"""
struct ModelPlatform <: Platform
    model   ::  Dictionary
end

model(p::ModelPlatform) = p.model
unsafe_dictionary(p::ModelPlatform, param) = resize(model(p), param)

param_first(p::ModelPlatform) = dimensions(model(p))
DictionaryStyle(p::ModelPlatform) = isbasis(model(p)) ? BasisStyle() : UnknownDictionaryStyle()
SamplingStyle(p::ModelPlatform) = isbasis(model(p)) ? InterpolationStyle() : OversamplingStyle()
SolverStyle(p::ModelPlatform, samplingstyle::SamplingStyle) = SolverStyle(DictionaryStyle(p), p, samplingstyle)

export correctparamformat
correctparamformat(p::ModelPlatform, param)  =
    param isa typeof(dimensions(model(p)))

export measure
measure(platform::ModelPlatform) = measure(model(platform))
coefficienttype(platform::ModelPlatform) = coefficienttype(model(platform))

components(platform::ModelPlatform) = map(ModelPlatform, components(model(platform)))
component(platform::ModelPlatform, i) = ModelPlatform(component(model(platform), i))

Base.complex(platform::ModelPlatform) = ModelPlatform(complex(platform.model))

#################
# ProductPlatform
#################

export ProductPlatform
"""
    struct ProductPlatform{N} <: Platform

A `ProductPlatform` corresponds to the product of two or more platforms. It results in a product
dictionary, product sampling operator, product solver, ... etcetera.
"""
struct ProductPlatform{N} <: Platform
    platforms   ::  NTuple{N,Platform}
end

ProductPlatform(platforms::Platform...) = ProductPlatform(platforms)

param_first(platform::ProductPlatform) = map(param_first, components(platform))

function DictionaryStyle(p::ProductPlatform)
    dict = TensorProductDict(map(x->dictionary(x,1), components(p))...)
    isbasis(dict) ? BasisStyle() : FrameStyle()
end
ProductPlatform(platform::Platform, n::Int) = ProductPlatform(ntuple(x->platform, n)...)
function dualdictionary(platform::ProductPlatform, param, measure::Measure; options...)
    @assert length(param)==length(components(platform))
    TensorProductDict(map((plati, parami, mi)->dualdictionary(plati, parami, mi; options...), components(platform), param, components(measure))...)
end
components(p::ProductPlatform) = p.platforms
factors(p::ProductPlatform) = components(p)

export productparameter
"""
    productparameter(p::ProductPlatform{N}, param)

Transform the parameter to a parameter accepted by a ProductPlatform
"""
productparameter(p::ProductPlatform{N}, n::Int) where {N} = ntuple(x->n, Val(N))
productparameter(p::ProductPlatform{N}, n::NTuple{N,Any}) where {N} = n

SamplingStyle(platform::ProductPlatform) = ProductSamplingStyle(map(SamplingStyle, components(platform)))
function SolverStyle(p::ProductPlatform, samplingstyle::ProductSamplingStyle)
    @assert nfactors(p) == ncomponents(samplingstyle)
    ProductSolverStyle(map(SolverStyle, factors(p), factors(samplingstyle)))
end

dictionary(p::ProductPlatform, n::Int) = dictionary(p, productparameter(p, n))

unsafe_dictionary(p::ProductPlatform, n) = tensorproduct(map(dictionary, components(p), n)...)

dualdictionary(platform::ProductPlatform, param, measure::Union{ProductWeight,DiscreteProductWeight}; options...) =
    tensorproduct(map((platformi, parami, mi)->dualdictionary(platformi, parami, mi; options...),
        components(platform), param, components(measure))...)

measure(platform::ProductPlatform) = productmeasure(map(measure, components(platform))...)

correctparamformat(p::ProductPlatform{N}, plt_par::NTuple{N,Int}) where N =
    all(map(correctparamformat, p.platforms, plt_par))

# we allow an integer platform parameter for convenience, but in that case
# the `productparameter` function gives the correct tuple
correctparamformat(p::ProductPlatform, plt_par::Int) = true

correctparamformat(::ProductPlatform, plt_par) = false

Base.complex(platform::ProductPlatform) =
    ProductPlatform(map(complex, components(platform)))


## Default platforms

export platformparameter, platform
"""
    platformparameter(dict::Dictionary) = dimensions(dict)

Return the parameter that is given to a model-based platform to obtain this dictionary.
"""
platformparameter(dict::Dictionary) = dimensions(dict)

"""
    platform(dict::Dictionary)

Return a platform that generates dictionaries of the type of dict.
"""
platform(dict::Dictionary) = ModelPlatform(dict)
platform(dict::TensorProductDict) = ProductPlatform(map(platform, factors(dict))...)
