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

import BasisFunctions: ComplexifiedDict

azdual_dict(dict::ComplexifiedDict; options...) =
    ComplexifiedDict(azdual_dict(superdict(dict); options...))
azdual_dict(dict::ExtensionFrame; options...) =
    ExtensionFrame(support(dict),azdual_dict(basis(dict); options...))

function azdual_dict(samplingstyle::ProductSamplingStyle, ap::ApproximationProblem; options...)
    TensorProductDict(   map( (x,style) ->azdual_dict(style, x; options...), elements(ap), samplingstyle.styles         )...      )
end

azdual_dict(samplingstyle::ProductSamplingStyle, ap::ApproximationProblem, measure::ProductMeasure; options...) =
    TensorProductDict(   map( (x,style, m) ->azdual_dict(style, x, m; options...), elements(ap), samplingstyle.styles, elements(measure)         )...      )
## Complexified dictionaries
