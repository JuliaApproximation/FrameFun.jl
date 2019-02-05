
# a dualplatformdictionary does not always need a measure (depending on the discretizationstyle), a dualdictionary does.
"The dual which is used to create a AZ `Z` matrix."
dualplatformdictionary(ap::ApproximationProblem; samplingstyle = SamplingStyle(ap), options...) =
    dualplatformdictionary(samplingstyle, ap; options...)

dualplatformdictionary(::DiscreteStyle, ap::DictionaryApproximation; options...) =
    dualplatformdictionary(dictionary(ap); options...)

dualplatformdictionary(::DiscreteStyle, ap::PlatformApproximation; options...) =
    dualplatformdictionary(ap.platform, ap.param; options... )

dualplatformdictionary(platform::Platform, n; dict = dictionary(platform, n), options...) =
    dualplatformdictionary(dict; options...)


function dualplatformdictionary(::SamplingStyle, ap::ApproximationProblem; options...)
    if iskey(options, :measure)
        measure = options[:measure]
    else
        measure = measure(ap; options...)
    end
    dualplatformdictionary(ap, measure; options...)
end
dualplatformdictionary(ap::DictionaryApproximation, measure::Measure; options...) =
    dualdictionary(ap, measure; options...)

dualplatformdictionary(ap::PlatformApproximation, measure::Measure; options...) =
    dualplatformdictionary(ap.platform, ap.param, measure; options... )

dualplatformdictionary(platform::Platform, n, measure::Measure; dict=platform(dict, n), options...) =
    dualdictionary(dict, measure; options...)

measure(ap::ApproximationProblem; options...) =
    haskey(options,:measure) ? options[:measure] : measure(ap)
measure(ap::DictionaryApproximation; options...) =
    haskey(options,:measure) ? options[:measure] : measure(dictionary(ap); options...)
measure(ap::PlatformApproximation; options...) =
    haskey(options,:measure) ? options[:measure] : measure(platform(ap), ap.param; options...)
measure(platform::Platform, n; dict = dictionary(platform, n), options...) =
    measure(dict; options...)







discrete_normalization(dict::Dictionary, L; S) = quadraturenormalization(S, measure(dict))

dualplatformdictionary(dict::Dictionary; options...) = _dualdictionary(dict, gramoperator(dict))
_dualdictionary(dict::Dictionary, gram::IdentityOperator) = dict
_dualdictionary(dict::Dictionary, gram) = conj(inv(gram)) * dict


## Weighted dictionaries

function dualplatformdictionary(dict::WeightedDict)
    dtilde = dualplatformdictionary(superdict(dict))
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

dualplatformdictionary(dict::TensorProductDict) = TensorProductDict(map(dualplatformdictionary, elements(dict))...)


## Complexified dictionaries

import BasisFunctions: ComplexifiedDict

dualplatformdictionary(dict::ComplexifiedDict) =
    ComplexifiedDict(dualplatformdictionary(superdict(dict)))


## Extension frames

dualplatformdictionary(dict::ExtensionFrame) = ExtensionFrame(support(dict),dualplatformdictionary(basis(dict)))
