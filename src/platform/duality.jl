
discrete_normalization(dict::Dictionary, L; measure=measure(dict), S) = quadraturenormalization(S, measure)

dualdictionary(dict::Dictionary) = _dualdictionary(dict, gramoperator(dict))
_dualdictionary(dict::Dictionary, gram::IdentityOperator) = dict
_dualdictionary(dict::Dictionary, gram) = conj(inv(gram)) * dict


## Weighted dictionaries

function dualdictionary(dict::WeightedDict)
    dtilde = dualdictionary(superdict(dict))
    invweight = x -> conj(1/dict.weightfun(x))
    WeightedDict(dtilde, invweight)
end

## Operated dictionaries

# function dualdictionary(dict::OperatedDict)
#     G = gramoperator(superdict(dict))
#     op = operator(dict)
#     OperatedDict(conj(inv(operator(dict))))
# end


## Mapped dictionaries

# dualdictionary(dict::MappedDict) = MappedDict(dualdictionary(superdict(dict)), mapping(dict))


## Product dictionaries

dualdictionary(dict::TensorProductDict) = TensorProductDict(map(dualdictionary, elements(dict))...)


## Complexified dictionaries

import BasisFunctions: ComplexifiedDict

dualdictionary(dict::ComplexifiedDict) =
    ComplexifiedDict(dualdictionary(superdict(dict)))


## Extension frames

dualdictionary(dict::ExtensionFrame) = ExtensionFrame(support(dict),dualdictionary(basis(dict)))
