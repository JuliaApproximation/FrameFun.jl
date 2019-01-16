
"""
A platform represents a family of dictionaries.

A platform typically has a primal and a dual sequence. The platform maps an
index to a set of parameter values, which is then used to generate a primal
or dual dictionary.
"""
abstract type Platform end

abstract type BasisPlatform <: Platform end
abstract type FramePlatform <: Platform end


"A dictionary style trait for platforms"
abstract type DictionaryStyle end

struct UnknownDictionaryStyle <: DictionaryStyle end
struct BasisStyle <: DictionaryStyle end
struct FrameStyle <: DictionaryStyle end

DictionaryStyle(p::Platform) = UnknownDictionaryStyle()
DictionaryStyle(p::BasisPlatform) = BasisStyle()
DictionaryStyle(p::FramePlatform) = FrameStyle()


"A sampling style trait"
abstract type SamplingStyle end

abstract type DiscreteStyle <: SamplingStyle end

"The sampling grid is determined by invoking `interpolation_grid` on the dictionary or platform."
struct InterpolationStyle <: DiscreteStyle end
"The sampling grid is determined by invoking `oversampling_grid` on the dictionary or platform."
struct OversamplingStyle <: DiscreteStyle end
"The sampling grid is determined by invoking `platform_grid` on the platform."
struct GridStyle <: DiscreteStyle end

abstract type ProjectionStyle <: SamplingStyle end
struct GramStyle <: ProjectionStyle end
struct RectangularGramStyle <: ProjectionStyle end

struct GenericSamplingStyle <: SamplingStyle end
struct ProductSamplingStyle <: SamplingStyle end


"A solver style trait"
abstract type SolverStyle end

struct DirectStyle <: SolverStyle end
struct InverseStyle <: SolverStyle end
struct DualStyle <: SolverStyle end
struct IterativeStyle <: SolverStyle end
struct TridiagonalProlateStyle <: SolverStyle end
struct AZStyle <: SolverStyle end
struct AZSmoothStyle <: SolverStyle end
struct ProductSolverStyle <: SolverStyle end


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
SamplingStyle(dict::TensorProductDict) = ProductSamplingStyle()


SolverStyle(dict::Dictionary, dstyle::SamplingStyle) = DirectStyle()
SolverStyle(dict::Dictionary, dstyle::InterpolationStyle) = hastransform(dict) ? InverseStyle() : DirectStyle()
SolverStyle(dict::TensorProductDict, dstyle::SamplingStyle) = ProductSolverStyle()
SolverStyle(dict::TensorProductDict, dstyle::InterpolationStyle) = ProductSolverStyle()
