module AugmentationPlatforms
using ..Platforms
using  BasisFunctions: Dictionary, Measure, MultiDict
import ..Platforms: SolverStyle, dictionary, dualdictionary, measure, param_next,
    correctparamformat, unsafe_dictionary

export AugmentationPlatform
"""
    struct AugmentationPlatform <: FramePlatform

An augmentation frame is the union of a basis, augmented with `K` extra
functions.

This type of frame is used as an example in the paper
"Frames and numerical approximation II: generalized sampling"
by B. Adcock and D. Huybrechs.
"""
struct AugmentationPlatform <: FramePlatform
    basis           ::  Platform
    functions       ::  Dictionary
end

correctparamformat(platform::AugmentationPlatform, param) =
    correctparamformat(platform.basis, param)

SolverStyle(platform::AugmentationPlatform, ::SamplingStyle) = AZStyle()

unsafe_dictionary(platform::AugmentationPlatform, i::Int) =
    MultiDict([unsafe_dictionary(platform.basis, i), platform.functions])
dualdictionary(platform::AugmentationPlatform, param, measure::Measure) =
    error("`dualdictionary` of `AugmentationPlatform` Not implemented.")
    # dualdictionary(platform.basis, param, measure)

measure(platform::AugmentationPlatform) = measure(platform.basis)

param_next(p::AugmentationPlatform, n::Int) = 2n

export ONB_plus_K
ONB_plus_K(basis::Platform, weightfunction, K) =
    AugmentationPlatform(basis, weightfunction * dictionary(basis, K))
end
