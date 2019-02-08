
# a dualplatformdictionary does not always need a measure (depending on the discretizationstyle), a dualdictionary does.
"""
The dual which is used to create a AZ `Z` matrix.

Can be overwritten at the level of the platform and the dictionary;
but it needs a specific implementation for `DiscreteStyle` since this
sampling style does not assume the implementation of `measure`
"""
@inline dualplatformdictionary(ap::ApproximationProblem; samplingstyle = SamplingStyle(ap), options...) =
    dualplatformdictionary(samplingstyle, ap; options...)

@inline dualplatformdictionary(samplingstyle::DiscreteStyle, ap::DictionaryApproximation; options...) =
    dualplatformdictionary(samplingstyle, dictionary(ap); options...)

@inline dualplatformdictionary(samplingstyle::DiscreteStyle, ap::PlatformApproximation; options...) =
    dualplatformdictionary(samplingstyle, ap.platform, ap.param; options...)

@inline dualplatformdictionary(samplingstyle::DiscreteStyle, platform::Platform, param;
            dict=dictionary(platform, param), options...) =
    dualplatformdictionary(samplingstyle, dict; dict=dict, options...)

@inline dualplatformdictionary(samplingstyle::SamplingStyle, ap::ApproximationProblem; options...) =
    dualplatformdictionary(samplingstyle, ap, measure(samplingstyle, ap; options...); options...)

@inline dualplatformdictionary(samplingstyle::SamplingStyle, ap::ApproximationProblem, measure::Measure; options...) =
    dualdictionary(dictionary(ap), measure; options...)

measure(ap::ApproximationProblem; samplingstyle, options...) =
    haskey(options,:measure) ? options[:measure] : measure(samplingstyle, ap; options...)
measure(samplingstyle::SamplingStyle, ap::ApproximationProblem; options...) =
    haskey(options,:measure) ? options[:measure] : measure(dictionary(ap); options...)
# Possibility on stackoverflow, but if neither measure with options or without options exists, this is an indication of an other problem.
# Other possibility. measure can never have optional arguments.
using InteractiveUtils
function measure(dict::Dictionary; options...)
    method = @which(measure(dict::Dictionary))
    if (String(method.file) == @__FILE__()) && (method.line == @__LINE__()-1)
        error("Implement `measure(::$(typeof(dict)))`")
    end
    measure(dict)
end

restrict(measure::BasisFunctions.DiscreteMeasure, domain::Domain) =
    BasisFunctions.DiscreteSubMeasure(subgrid(grid(measure), domain), weights(measure))

measure(::DiscreteGramStyle, ap::ApproximationProblem; options...) = discrete_gram_measure(ap; options...)


"""


Can be specialized for platforms and dictionaries a the level of `discrete_gram_measure`
or at the level of `oversampling_grid`, `oversampling_weight`
"""
function discrete_gram_measure(ap; options...)
    L = samplingparameter(OversamplingStyle(), ap; options...)
    discrete_gram_measure(ap, L; options...)
end

discrete_gram_measure(ap::DictionaryApproximation, L; options...) =
    discrete_gram_measure(dictionary(ap), L)
discrete_gram_measure(ap::PlatformApproximation, L; options...) =
    discrete_gram_measure(ap.platform, ap.param, L; options...)
discrete_gram_measure(platform::Platform, n, L=n; dict = dictionary(platform, n), options...) =
    discrete_gram_measure(dict, L)

function discrete_gram_measure(dict::Dictionary, L)
    grid = oversampling_grid(dict, L)
    weight = oversampling_weight(dict, grid)
    BasisFunctions.DiscreteMeasure(grid, weight)
end

oversampling_weight(dict::Dictionary, grid::AbstractGrid) = Ones{coefficienttype(dict)}(size(grid)...)
oversampling_weight(dict::Dictionary, grid::AbstractSubGrid) = Ones{coefficienttype(dict)}(size(supergrid(grid))...)

@inline dualplatformdictionary(::DiscreteStyle, dict::Dictionary; options...) =
    dualplatformdictionary(dict; options...)

@inline dualplatformdictionary(dict::Dictionary; options...) = _dualdictionary(dict, gramoperator(dict))
_dualdictionary(dict::Dictionary, gram::IdentityOperator) = dict
_dualdictionary(dict::Dictionary, gram) = conj(inv(gram)) * dict


## Weighted dictionaries

function dualplatformdictionary(dict::WeightedDict; options...)
    dtilde = dualplatformdictionary(superdict(dict); options...)
    invweight = x -> conj(1/dict.weightfun(x))
    WeightedDict(dtilde, invweight)
end

## Operated dictionaries

# function dualplatformdictionary(dict::OperatedDict)
#     G = gramoperator(superdict(dict))
#     op = operator(dict)
#     OperatedDict(conj(inv(operator(dict))))
# end


## Mapped dictionaries

# dualplatformdictionary(dict::MappedDict) = MappedDict(dualplatformdictionary(superdict(dict)), mapping(dict))


## Product dictionaries

dualplatformdictionary(dict::TensorProductDict; options...) = TensorProductDict(map(dualplatformdictionary, elements(dict))...)


## Complexified dictionaries

import BasisFunctions: ComplexifiedDict

dualplatformdictionary(dict::ComplexifiedDict; options...) =
    ComplexifiedDict(dualplatformdictionary(superdict(dict); options...))


## Extension frames

dualplatformdictionary(dict::ExtensionFrame; options...) = ExtensionFrame(support(dict),dualplatformdictionary(basis(dict); options...))
