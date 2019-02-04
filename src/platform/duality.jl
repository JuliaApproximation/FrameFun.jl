
discrete_normalization(dict::Dictionary, L; measure=measure(dict), S) = quadraturenormalization(S, measure)


## Complexified dictionaries

import BasisFunctions: ComplexifiedDict

dualdictionary(dict::ComplexifiedDict, measure::BasisFunctions.Measure=measure(dict); options...) =
    ComplexifiedDict(dualdictionary(superdict(dict), measure; options...))


## Extension frames

dualdictionary(dict::ExtensionFrame, measure::BasisFunctions.Measure=measure(dict)) =
    ExtensionFrame(support(dict),dualdictionary(basis(dict), BasisFunctions.supermeasure(measure)))
