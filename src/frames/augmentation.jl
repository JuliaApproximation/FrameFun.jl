
"""
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

weightfun(platform::AugmentationPlatform) = platform.weightfun

dictionary(platform::AugmentationPlatform, i::Int) =
    MultiDict([dictionary(platform.basis, i), platform.functions])

dualplatformdictionary(sstyle::DiscreteStyle, platform::AugmentationPlatform, i::Int; options...) =
    MultiDict([dualplatformdictionary(sstyle, platform.basis, i; options...), platform.functions])

oversampling_grid(platform::AugmentationPlatform, param, L; dict, options...) =
    oversampling_grid(platform.basis, param, L)

discrete_normalization(platform::AugmentationPlatform, n, L; options...) =
    discrete_normalization(platform.basis, n, L; options...)

measure(platform::AugmentationPlatform) = measure(platform.basis)

param_next(p::AugmentationPlatform, n::Int) = 2n

ONB_plus_K(basis::Platform, weightfunction, K) =
    AugmentationPlatform(basis, weightfunction * dictionary(basis, K))
