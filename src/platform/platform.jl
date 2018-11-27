
"""
A platform represents a family of dictionaries.

A platform typically has a primal and a dual sequence. The platform maps an
index to a set of parameter values, which is then used to generate a primal
or dual dictionary.
"""
abstract type Platform end

abstract type BasisPlatform <: Platform end
abstract type FramePlatform <: Platform end


"A dictionary style trait"
abstract type DictionaryStyle end

struct UnknownDictionaryStyle <: DictionaryStyle end
struct BasisStyle <: DictionaryStyle end
struct FrameStyle <: DictionaryStyle end

DictionaryStyle(p::Platform) = UnknownDictionaryStyle()
DictionaryStyle(p::BasisPlatform) = BasisStyle()
DictionaryStyle(p::FramePlatform) = FrameStyle()


"A solver style trait"
abstract type SolverStyle end

struct DirectStyle <: SolverStyle end
struct InverseStyle <: SolverStyle end
struct TransformStyle <: SolverStyle end
struct IterativeStyle <: SolverStyle end
struct AZStyle <: SolverStyle end
struct TridiagonalProlateStyle <: SolverStyle end

SolverStyle(p::Platform) = SolverStyle(DictionaryStyle(p), p)

SolverStyle(::BasisStyle, p::Platform) = InverseStyle()
SolverStyle(::FrameStyle, p::Platform) = AZStyle()
SolverStyle(::DictionaryStyle, p::Platform) = DirectStyle()


"A discretization style trait"
abstract type DiscretizationStyle end

abstract type DiscreteSamplingStyle <: DiscretizationStyle end
struct InterpolationStyle <: DiscreteSamplingStyle end
struct OversamplingStyle <: DiscreteSamplingStyle end
struct GridStyle <: DiscreteSamplingStyle end

abstract type ProjectionStyle <: DiscretizationStyle end
struct GramStyle <: ProjectionStyle end
struct RectangularGramStyle <: ProjectionStyle end

struct GenericSamplingStyle <: DiscretizationStyle end

DiscretizationStyle(p::Platform) = DiscretizationStyle(DictionaryStyle(p), p)

DiscretizationStyle(::BasisStyle, p::Platform) = InterpolationStyle()
DiscretizationStyle(::FrameStyle, p::Platform) = OverSamplingStyle()
