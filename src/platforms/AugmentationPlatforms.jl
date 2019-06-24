module AugmentationPlatforms
using ..Platforms
using  BasisFunctions: Dictionary, Measure, MultiDict
import ..Platforms: SolverStyle, dictionary, dualdictionary, measure, param_next

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

SolverStyle(platform::AugmentationPlatform, ::SamplingStyle) = AZStyle()

dictionary(platform::AugmentationPlatform, i::Int) =
    MultiDict([dictionary(platform.basis, i), platform.functions])
dualdictionary(platform::AugmentationPlatform, param, ::Measure) =
    error("Not implemented")

measure(platform::AugmentationPlatform) = measure(platform.basis)

param_next(p::AugmentationPlatform, n::Int) = 2n

export ONB_plus_K
ONB_plus_K(basis::Platform, weightfunction, K) =
    AugmentationPlatform(basis, weightfunction * dictionary(basis, K))
end
