
using  BasisFunctions: Dictionary, Measure, MultiDict

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
function dualdictionary(platform::AugmentationPlatform, param, measure::Measure; options...)
#    error("`dualdictionary` of `AugmentationPlatform` Not implemented.")
    # TODO: line below is a hack. Problem is you want to generate a dual of
    # the same size as the dictionary itself, but we only know the dual of its
    # basis component. We add the length of the extra functions and hope for the best.
    @warn "Routine 'dualdictionary' invoked on AugmentationPlatform does not know how to generically produce a dual of the right length."
    dualdictionary(platform.basis, param+length(platform.functions), measure)
end

measure(platform::AugmentationPlatform) = measure(platform.basis)

param_double(p::AugmentationPlatform, n::Int) = 2n

export ONB_plus_K
ONB_plus_K(basis::Platform, weightfunction, K) =
    AugmentationPlatform(basis, weightfunction * dictionary(basis, K))
