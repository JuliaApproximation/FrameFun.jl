"""
    dualplatformdictionary(ap::ApproximationProblem; options...)

The dual that is used to create a AZ `Z` matrix.
"""
dualplatformdictionary(ap::ApproximationProblem; samplingstyle = SamplingStyle(ap), options...) =
    dualplatformdictionary(samplingstyle, ap; options...)

dualplatformdictionary(samplingstyle::DiscreteStyle, ap::ApproximationProblem; options...) =
    dualplatformdictionary(samplingstyle, ap, discretemeasure(samplingstyle, ap); options...)

dualplatformdictionary(samplingstyle::SamplingStyle, ap::ApproximationProblem; options...) =
    dualplatformdictionary(samplingstyle, ap, measure(samplingstyle, ap; options...); options...)

dualplatformdictionary(samplingstyle::SamplingStyle, ap::ApproximationProblem, measure::Measure; options...) =
    azdual(dictionary(ap), measure; options...)

azdual(dict::Dictionary, measure; options...) =
    BasisFunctions.default_gramdual(dict, measure; options...)

azdual(dict::ExtensionFrame, measure; options...) =
    extensionframe(support(dict), azdual(superdict(dict), supermeasure(measure); options...),)

function azdual(dict::MultiDict, measure::Measure; options...)
    dictionary = dict.dicts[1].superdict
    weights = map(weightfunction, elements(dict))
    denom = (x...)->sum(map(w->abs(w(x...))^2, weights))
    MultiDict([((x...)->(weights[j](x...)/denom(x...))) * azdual(dictionary, measure; options...) for j=1:length(weights)])
end

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

measure(ss::DiscreteGramStyle, ap::ApproximationProblem; options...) =
    discrete_gram_measure(ss, ap; options...)

"""
    @aptoplatform `function`

Implement the function at the level of the dictionary.
"""
macro aptoplatform(ex)
    def = Meta.parse("default_"*string(ex))
    ret = quote
        $(ex)(ss::SamplingStyle, ap::ApproximationProblem; options...) =
            $(ex)(ss, ap, samplingparameter(ss, ap; options...); options...)
        $(ex)(ss::SamplingStyle, ap::DictionaryApproximation, L; options...) =
            $(ex)(ss, dictionary(ap), L)
        $(ex)(ss::SamplingStyle, ap::PlatformApproximation, L; options...) =
            $(ex)(ss, ap.platform, ap.param, L; options...)
        $(ex)(ss::SamplingStyle, platform::Platform, n, L=n; dict = dictionary(platform, n), options...) =
            $(ex)(ss, dict, L)
        $(ex)(ss::SamplingStyle, dict::Dictionary, L; options...) =
            $(def)(ss, dict, L)
    end
    esc(ret)
end

@aptoplatform discrete_gram_measure

default_discrete_gram_measure(ss::DiscreteGramStyle, dict::Dictionary, L) =
    discretemeasure(oversampling_grid(dict, L))

@aptoplatform discretemeasure

default_discretemeasure(ss::DiscreteStyle, dict::Dictionary, L) =
    discretemeasure(oversampling_grid(dict, L))



dualplatformdictionary(dict::Dictionary; measure = measure(dict), options...) =
    _dual(dict, gramoperator(dict, measure))
_dual(dict::Dictionary, gram::IdentityOperator) = dict
_dual(dict::Dictionary, gram) = conj(inv(gram)) * dict


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

dualplatformdictionary(dict::TensorProductDict; options...) = TensorProductDict(map(x->dualplatformdictionary(x;options...), elements(dict))...)
function dualplatformdictionary(samplingstyle::ProductSamplingStyle, ap::ApproximationProblem; options...)
    TensorProductDict(   map( (x,style) ->dualplatformdictionary(style, x; options...), elements(ap), samplingstyle.styles         )...      )
end

dualplatformdictionary(samplingstyle::ProductSamplingStyle, ap::ApproximationProblem, measure::ProductMeasure; options...) =
    TensorProductDict(   map( (x,style, m) ->dualplatformdictionary(style, x, m; options...), elements(ap), samplingstyle.styles, elements(measure)         )...      )
## Complexified dictionaries

import BasisFunctions: ComplexifiedDict

dualplatformdictionary(dict::ComplexifiedDict; options...) =
    ComplexifiedDict(dualplatformdictionary(superdict(dict); options...))
dualplatformdictionary(dict::ExtensionFrame; options...) =
    ExtensionFrame(support(dict),dualplatformdictionary(basis(dict); options...))
