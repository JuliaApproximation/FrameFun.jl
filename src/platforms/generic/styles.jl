
#################
# DictionaryStyle
#################

export DictionaryStyle, UnknownDictionaryStyle, BasisStyle, FrameStyle
"""
    abstract type DictionaryStyle end

A dictionary style trait for platforms

See also, [`BasisStyle`](@ref), [`FrameStyle`](@ref),
[`UnknownDictionaryStyle`](@ref)
"""
abstract type DictionaryStyle end

"""
    struct UnknownDictionaryStyle <: DictionaryStyle end

The platform creates dictionaries of unknown format.
"""
struct UnknownDictionaryStyle <: DictionaryStyle end

"""
    struct BasisStyle <: DictionaryStyle end

The platform creates dictionaries that are bases.
"""
struct BasisStyle <: DictionaryStyle end

"""
    struct FrameStyle <: DictionaryStyle end

The platform creates dictionaries that are frames.
"""
struct FrameStyle <: DictionaryStyle end

DictionaryStyle(platform::Platform) = UnknownDictionaryStyle()
DictionaryStyle(platform::BasisPlatform) = BasisStyle()
DictionaryStyle(platform::FramePlatform) = FrameStyle()







#################
# SamplingStyle
#################

export SamplingStyle, DiscreteStyle, ProjectionStyle, InterpolationStyle,
    OversamplingStyle, GridStyle, GramStyle, DiscreteGramStyle,
    RectangularGramStyle, GenericSamplingStyle, ProductSamplingStyle
"""
    abstract type SamplingStyle end

A sampling style trait for platforms

See also [`DiscreteStyle`](@ref), [`ProjectionStyle`](@ref),
[`GenericSamplingStyle`](@ref), [`ProductSamplingStyle`](@ref)
"""
abstract type SamplingStyle end

"""
    abstract type DiscreteStyle <: SamplingStyle end

The sampling operator corresponds to evaluation in a grid. This grid is
specified by the listed styles here:

[`InterpolationStyle`](@ref), [`OversamplingStyle`](@ref),
[`GridStyle`](@ref)
"""
abstract type DiscreteStyle <: SamplingStyle end

"Discrete sampling without associated weights."
abstract type UnweightedDiscreteStyle <: DiscreteStyle end
"Discrete sampling with associated weights."
abstract type WeightedDiscreteStyle <: DiscreteStyle end

"The sampling grid is determined by invoking `interpolation_grid` on the dictionary or platform."
struct InterpolationStyle <: UnweightedDiscreteStyle end
"The sampling grid is determined by invoking `oversampling_grid` on the dictionary or platform."
struct OversamplingStyle <: UnweightedDiscreteStyle end
"""
The sampling grid is determined by invoking `platform_grid` on the platform,
or by providing a `grid=...` argument to the Fun constructor.
"""
struct GridStyle <: UnweightedDiscreteStyle end

"Discrete sampling followed by weighting the samples."
struct WeightedSampling <: WeightedDiscreteStyle
    style       # the discrete sampling style
end
unweighted(st::WeightedSampling) = st.style

"Discrete sampling followed by normalizing the samples."
struct NormalizedSampling <: WeightedDiscreteStyle
    style       # the discrete sampling style
end
unweighted(st::NormalizedSampling) = st.style

NormalizedSampling(style::NormalizedSampling) = style


"""
    abstract type ProjectionStyle <: SamplingStyle end

The sampling operator corresponds to projections (using inner products).
The kind of inner product is specified by
[`GramStyle`](@ref), [`DiscreteGramStyle`](@ref),
[`RectangularGramStyle`](@ref)
"""
abstract type ProjectionStyle <: SamplingStyle end
"""
    struct GramStyle <: ProjectionStyle end

Use the inner product implied by the measure obtained by invoking `measure` on the
dictionary or platform. This results in a Gram matrix.
"""
struct GramStyle <: ProjectionStyle end
"""
    struct DiscreteGramStyle <: ProjectionStyle end

Use the inner product implied by the measure obtained by invoking `discretemeasure` on the
dictionary or platform. This results in a discrete Gram matrix.
"""
struct DiscreteGramStyle <: ProjectionStyle end
"""
    struct RectangularGramStyle <: ProjectionStyle end

Not implemented (oversample with a projection on a larger dual dictionary)
"""
struct RectangularGramStyle <: ProjectionStyle end

"""
    struct GenericSamplingStyle <: SamplingStyle end

Write your own sampling_operator. The sampling operator corresponds the one implemented
by the user.
"""
struct GenericSamplingStyle <: SamplingStyle end

"""
    struct ProductSamplingStyle{STYLES} <: SamplingStyle

The sampling operator has product structure.
"""
struct ProductSamplingStyle{STYLES} <: SamplingStyle
    styles  :: STYLES
end
ProductSamplingStyle(styles::SamplingStyle...) = ProductSamplingStyle(styles)
components(style::ProductSamplingStyle) = style.styles
factors(style::ProductSamplingStyle) = components(style)

# TODO: remove this definition
export SamplingStyleSuper
const SamplingStyleSuper{TYPE} = Union{<:TYPE,ProductSamplingStyle{NTuple{N,<:TYPE}} where N} where {TYPE<:SamplingStyle}




#################
# SolverStyle
#################

export SolverStyle, DirectStyle, InverseStyle, DualStyle, IterativeStyle,
    TridiagonalProlateStyle, AZStyle, AZSmoothStyle, ProductSolverStyle
"""
    abstract type SolverStyle end

A solver style trait for platforms

See also [`DirectStyle`](@ref), [`InverseStyle`](@ref),
[`DualStyle`](@ref), [`IterativeStyle`](@ref),
[`TridiagonalProlateStyle`](@ref), [`AZStyle`](@ref),
[`AZSmoothStyle`](@ref), [`ProductSolverStyle`](@ref)
"""
abstract type SolverStyle end
"""
    struct DirectStyle <: SolverStyle end

Use a direct solver to solve the discretized approximation problem.
"""
struct DirectStyle <: SolverStyle end
"""
    struct InverseStyle <: SolverStyle end

Use the inverse operator to solve the discretized approximation problem.
"""
struct InverseStyle <: SolverStyle end
"""
    struct DualStyle <: SolverStyle end

Use projection on the dual dictionary to solve the discretized approximation problem.
"""
struct DualStyle <: SolverStyle end
"""
    struct IterativeStyle <: SolverStyle end

Use an iterative method to solve the discretized approximation problem.
"""
struct IterativeStyle <: SolverStyle end
"""
    struct TridiagonalProlateStyle <: SolverStyle end

Use a tridiagonal solver to solve the discretized approximation problem
(only for Fourier extension)
"""
struct TridiagonalProlateStyle <: SolverStyle end
"""
    struct AZStyle <: SolverStyle end

Use the AZ algorithm to solve the discretized approximation problem.
"""
struct AZStyle <: SolverStyle end
"""
    struct AZSmoothStyle <: SolverStyle end

Use the AZ algorithm to solve the discretized approximation problem,
but tweak the coefficients such that a smooth extension is obtained.
"""
struct AZSmoothStyle <: SolverStyle end

"""
    struct ProductSolverStyle <: SolverStyle

Use a tensor product solver to solve the discretized approximation problem.
"""
struct ProductSolverStyle <: SolverStyle
    styles
end
ProductSolverStyle(styles::SolverStyle...) = ProductSolverStyle(styles)
components(style::ProductSolverStyle) = style.styles
factors(style::ProductSolverStyle) = components(style)



#################
# Defaults
#################

# Defaults for platforms
SamplingStyle(p::Platform) = SamplingStyle(DictionaryStyle(p), p)

SamplingStyle(::BasisStyle, p::Platform) = InterpolationStyle()
SamplingStyle(::FrameStyle, p::Platform) = OversamplingStyle()
SamplingStyle(::UnknownDictionaryStyle, p::Platform) = OversamplingStyle()


SolverStyle(p::Platform, samplingstyle) =
    SolverStyle(DictionaryStyle(p), p, samplingstyle)

SolverStyle(::BasisStyle, p::Platform, ::InterpolationStyle) = InverseStyle()
SolverStyle(::BasisStyle, p::Platform, ::OversamplingStyle) = DirectStyle()
SolverStyle(::BasisStyle, p::Platform, ::GridStyle) = DirectStyle()
SolverStyle(::BasisStyle, p::Platform, ::SamplingStyle) = DirectStyle()
SolverStyle(::FrameStyle, p::Platform, ::SamplingStyle) = AZStyle()
SolverStyle(::UnknownDictionaryStyle, p::Platform, ::SamplingStyle) = DirectStyle()
