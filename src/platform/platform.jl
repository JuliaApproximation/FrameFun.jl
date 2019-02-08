
"""
A platform represents a family of dictionaries.

A platform typically has a primal and a dual sequence. The platform maps an
index to a set of parameter values, which is then used to generate a primal
or dual dictionary.
"""
abstract type Platform end

"""
A `BasisPlatform` represents a family of bases.
"""
abstract type BasisPlatform <: Platform end

"""
A `FramePlatform` represents a family of frames
(or ill-conditioned bases which are in the limit a frame).
"""
abstract type FramePlatform <: Platform end


"A dictionary style trait for platforms"
abstract type DictionaryStyle end

struct UnknownDictionaryStyle <: DictionaryStyle end
struct BasisStyle <: DictionaryStyle end
struct FrameStyle <: DictionaryStyle end

DictionaryStyle(platform::Platform) = UnknownDictionaryStyle()
DictionaryStyle(platform::BasisPlatform) = BasisStyle()
DictionaryStyle(platform::FramePlatform) = FrameStyle()

getindex(platform::Platform, n) = dictionary(platform, n)

"A sampling style trait"
abstract type SamplingStyle end

"The sampling operator corresponds to evaluation in a discrete set of points."
abstract type DiscreteStyle <: SamplingStyle end

"The sampling grid is determined by invoking `interpolation_grid` on the dictionary or platform."
struct InterpolationStyle <: DiscreteStyle end
"The sampling grid is determined by invoking `oversampling_grid` on the dictionary or platform."
struct OversamplingStyle <: DiscreteStyle end
"""
The sampling grid is determined by invoking `platform_grid` on the platform, or
by providing a `grid=...` argument to the Fun constructor.
"""
struct GridStyle <: DiscreteStyle end


"The sampling operator corresponds to projections (using inner products)."
abstract type ProjectionStyle <: SamplingStyle end
struct GramStyle <: ProjectionStyle end
struct DiscreteGramStyle <: ProjectionStyle end
struct RectangularGramStyle <: ProjectionStyle end

"The sampling operator is determined by invoking `genericsamplingoperator` on the platform."
struct GenericSamplingStyle <: SamplingStyle end

struct ProductSamplingStyle <: SamplingStyle
    styles
end
ProductSamplingStyle(styles::SamplingStyle...) = ProductSamplingStyle(styles)


"A solver style trait"
abstract type SolverStyle end

struct DirectStyle <: SolverStyle end
struct InverseStyle <: SolverStyle end
struct DualStyle <: SolverStyle end
struct IterativeStyle <: SolverStyle end
struct TridiagonalProlateStyle <: SolverStyle end
struct AZStyle <: SolverStyle end
struct AZSmoothStyle <: SolverStyle end

struct ProductSolverStyle <: SolverStyle
    styles
end
ProductSolverStyle(styles::SolverStyle...) = ProductSolverStyle(styles)

## Defaults

# Defaults for platforms
SamplingStyle(p::Platform) = SamplingStyle(DictionaryStyle(p), p)

SamplingStyle(::BasisStyle, p::Platform) = InterpolationStyle()
SamplingStyle(::FrameStyle, p::Platform) = OversamplingStyle()


SolverStyle(p::Platform, dstyle::SamplingStyle) =
    SolverStyle(DictionaryStyle(p), p, dstyle)

SolverStyle(::BasisStyle, p::Platform, ::InterpolationStyle) = InverseStyle()
SolverStyle(::BasisStyle, p::Platform, ::OversamplingStyle) = DirectStyle()
SolverStyle(::BasisStyle, p::Platform, ::GridStyle) = DirectStyle()
SolverStyle(::BasisStyle, p::Platform, dstyle) = DirectStyle()
SolverStyle(::FrameStyle, p::Platform, dstyle) = AZStyle()
SolverStyle(::UnknownDictionaryStyle, p::Platform, dstyle) = DirectStyle()


# Defaults for dictionaries
SamplingStyle(dict::Dictionary) = isbasis(dict) ? InterpolationStyle() : OversamplingStyle()
SamplingStyle(dict::TensorProductDict) = ProductSamplingStyle(map(SamplingStyle, elements(dict)))


SolverStyle(dict::Dictionary, dstyle::SamplingStyle) = DirectStyle()
SolverStyle(dict::Dictionary, dstyle::InterpolationStyle) = hastransform(dict) ? InverseStyle() : DirectStyle()
SolverStyle(dict::TensorProductDict, dstyle::ProductSamplingStyle) =
    ProductSolverStyle(map(SolverStyle, elements(dict), dstyle.styles))
