
# Determine the measure
"""
    azdual_dict(ap::ApproximationProblem; options...)

The dual that is used to create a AZ `Z` matrix.
"""
azdual_dict(samplingstyle::DiscreteStyle, ap::ApproximationProblem; options...) =
    azdual_dict(samplingstyle, ap, discretemeasure(samplingstyle, ap; options...); options...)
azdual_dict(samplingstyle::SamplingStyle, ap::ApproximationProblem; options...) =
    azdual_dict(samplingstyle, ap, measure(samplingstyle, ap; options...); options...)

azdual_dict(ss::ProductSamplingStyle, ap::ApproximationProblem; options...) =
    TensorProductDict( ap_components(azdual_dict, ss, ap; options...)... )

azdual_dict(samplingstyle::SamplingStyle, ap::Platform, param, measure; options...) =
    dualdictionary(ap, param, measure; samplingstyle=samplingstyle, options...)
default_azdual_dict(::SamplingStyle, dict::Dictionary, measure; options...) =
    (@warn "azdualdict on dictionary is deprecated";BasisFunctions.default_gramdual(dict, measure; options...))
