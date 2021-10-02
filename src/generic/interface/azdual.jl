
import BasisFunctions: DualType

"Trait type to select a dual dictionary suitable for the AZ algorithm."
struct AZDual <: DualType end

dual(::AZDual, dict::Dictionary, measure; options...) =
    azdual(dict, measure; options...)

"""
The dual that is used to create an AZ `Z` matrix.
"""
azdual(ss::DiscreteStyle, ap::ApproximationProblem; options...) =
    dualdictionary(platform(ap), platformparameter(ap), discretemeasure(ap); samplingstyle=ss, options...)
azdual(ss::SamplingStyle, ap::ApproximationProblem; options...) =
    dualdictionary(platform(ap), platformparameter(ap), measure(ap); samplingstyle=ss, options...)
