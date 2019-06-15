"""
    @aptoplatform `function`

Implement the function at the level of the dictionary.
"""
macro aptoplatform(ex)
    def = Meta.parse("default_"*string(ex))
    intermediate = Meta.parse("_"*string(ex))
    ret = quote
        $(ex)(ap::ApproximationProblem, args...; samplingstyle=SamplingStyle(ap),options...) =
            $(intermediate)(samplingstyle, ap, samplingparameter(samplingstyle, ap; options...), args...; options...)
        $(ex)(ss::SamplingStyle, ap::ApproximationProblem, args...; options...) =
            $(intermediate)(ss, ap, samplingparameter(ss, ap; options...), args...; options...)
        $(intermediate)(ss::SamplingStyle, ap::DictionaryApproximation, L, args...; options...) =
            $(ex)(ss, dictionary(ap), L, args...)
        $(intermediate)(ss::SamplingStyle, ap::PlatformApproximation, L, args...; options...) =
            $(ex)(ss, ap.platform, ap.param, L, args...; options...)
        $(ex)(ss::SamplingStyle, platform::Platform, param, L, args...; dict = dictionary(platform, param), options...) =
            $(ex)(ss, dict, L, args...)
        $(ex)(ss::SamplingStyle, dict::Dictionary, L, args...; options...) =
            $(def)(ss, dict, L, args...)
    end
    esc(ret)
end

"""
    azdual_dict(ap::ApproximationProblem; options...)

The dual that is used to create a AZ `Z` matrix.
"""
azdual_dict(samplingstyle::DiscreteStyle, ap::ApproximationProblem; options...) =
    azdual_dict(samplingstyle, ap, discretemeasure(samplingstyle, ap); options...)

azdual_dict(samplingstyle::SamplingStyle, ap::ApproximationProblem; options...) =
    azdual_dict(samplingstyle, ap, measure(samplingstyle, ap; options...); options...)

azdual_dict(ap::ApproximationProblem; samplingstyle=SamplingStyle(ap), options...) =
    azdual_dict(samplingstyle, ap; options...)

@aptoplatform azdual_dict

default_azdual_dict(::SamplingStyle, dict::Dictionary, L, measure::Measure; options...) =
    BasisFunctions.default_gramdual(dict, measure; options...)
azdual_dict(sstyle::SamplingStyle, dict::ExtensionFrame, L, measure::Measure; options...) =
   extensionframe(support(dict), azdual_dict(sstyle, superdict(dict), L, supermeasure(measure); options...),)

function azdual_dict(dict::MultiDict, measure::Measure; options...)
    dictionary = dict.dicts[1].superdict
    weights = map(weightfunction, elements(dict))
    denom = (x...)->sum(map(w->abs(w(x...))^2, weights))
    MultiDict([((x...)->(weights[j](x...)/denom(x...))) * azdual_dict(dictionary, measure; options...) for j=1:length(weights)])
end

@aptoplatform measure
default_measure(::SamplingStyle, dict::Dictionary, L; options...) =
    measure(dict)

measure(ss::DiscreteGramStyle, ap::ApproximationProblem; options...) =
    discrete_gram_measure(ss, ap; options...)

@aptoplatform discrete_gram_measure

default_discrete_gram_measure(ss::DiscreteGramStyle, dict::Dictionary, L) =
    discretemeasure(oversampling_grid(dict, L))

@aptoplatform discretemeasure

default_discretemeasure(ss::DiscreteStyle, dict::Dictionary, L) =
    discretemeasure(oversampling_grid(dict, L))

azdual_dict(dict::Dictionary; measure = measure(dict), options...) =
    _dual(dict, gramoperator(dict, measure))
_dual(dict::Dictionary, gram::IdentityOperator) = dict
_dual(dict::Dictionary, gram) = conj(inv(gram)) * dict


## Weighted dictionaries

function azdual_dict(dict::WeightedDict; options...)
    dtilde = azdual_dict(superdict(dict); options...)
    invweight = x -> conj(1/dict.weightfun(x))
    WeightedDict(dtilde, invweight)
end

## Operated dictionaries

# function azdual_dict(dict::OperatedDict)
#     G = gramoperator(superdict(dict))
#     op = operator(dict)
#     OperatedDict(conj(inv(operator(dict))))
# end


## Mapped dictionaries

# azdual_dict(dict::MappedDict) = MappedDict(azdual_dict(superdict(dict)), mapping(dict))


## Product dictionaries

azdual_dict(dict::TensorProductDict; options...) = TensorProductDict(map(x->azdual_dict(x;options...), elements(dict))...)
function azdual_dict(samplingstyle::ProductSamplingStyle, ap::ApproximationProblem; options...)
    TensorProductDict(   map( (x,style) ->azdual_dict(style, x; options...), elements(ap), samplingstyle.styles         )...      )
end

azdual_dict(samplingstyle::ProductSamplingStyle, ap::ApproximationProblem, measure::ProductMeasure; options...) =
    TensorProductDict(   map( (x,style, m) ->azdual_dict(style, x, m; options...), elements(ap), samplingstyle.styles, elements(measure)         )...      )
## Complexified dictionaries

import BasisFunctions: ComplexifiedDict

azdual_dict(dict::ComplexifiedDict; options...) =
    ComplexifiedDict(azdual_dict(superdict(dict); options...))
azdual_dict(dict::ExtensionFrame; options...) =
    ExtensionFrame(support(dict),azdual_dict(basis(dict); options...))
