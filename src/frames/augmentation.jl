
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

dualdictionary(platform::AugmentationPlatform, i::Int; dict = dictionary(platform, i)) =
    MultiDict([dualdictionary(platform.basis, i), platform.functions])

oversampling_grid(platform::AugmentationPlatform, param, L; dict, options...) =
    oversampling_grid(platform.basis, param, L)

discrete_normalization(platform::AugmentationPlatform, n, L; options...) =
    discrete_normalization(platform.basis, n, L; options...)

measure(platform::AugmentationPlatform) = measure(platform.basis)

nextsize(p::AugmentationPlatform, n::Int) = 2n

ONB_plus_K(basis::Platform, weightfunction, K) =
    AugmentationPlatform(basis, weightfunction * dictionary(basis, K))
