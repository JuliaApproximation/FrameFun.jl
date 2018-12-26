
dualsamplingoperator(dict::Dictionary; options...) = dualsamplingoperator(dict, samplingoperator(dict; options...); options...)

dualsamplingoperator(dict::Dictionary, m::Int; options...) =
    dualsamplingoperator(dict, samplingoperator(dict, M=m; options...))

dualsamplingoperator(dict::Dictionary, S) =
    quadraturenormalization(S) * S


dualdictionary(dict::FourierBasis) = dict

function dualdictionary(dict::ChebyshevBasis{T}) where {T}
    scaling = ScalingOperator(dict, 2/convert(T, pi)) *
        BasisFunctions.CoefficientScalingOperator(dict, 1, one(T)/2)
    scaling * dict
end

## Weighted dictionaries

function dualdictionary(dict::WeightedDict)
    dtilde = dualdictionary(superdict(dict))
    invweight = x -> 1/dict.weightfun(x)
    WeightedDict(dtilde, invweight)
end

# dualsamplingoperator(dict::WeightedDictionary, m; T = coefficienttype(dict), S = samplingoperator(dict, m, T=T)) =
    # quadraturenormalization(S) * S

## Operated dictionaries

function dualdictionary(dict::OperatedDict)
    op = operator(dict)
    @assert is_diagonal(op)
    invop = similar(inv(op), dualdictionary(src(op)), dualdictionary(dest(op)))
    OperatedDict(invop)
end


## Mapped dictionaries

dualdictionary(dict::MappedDict) = MappedDict(dualdictionary(superdict(dict)), mapping(dict))

dualsamplingoperator(dict::MappedDict, S) =
    dualsamplingoperator(superdict(dict), S)


## Product dictionaries

dualdictionary(dict::TensorProductDict) = TensorProductDict(map(dualdictionary, elements(dict))...)

function dualsamplingoperator(dict::TensorProductDict, S::GridSampling)
    duals = map(dualsamplingoperator, elements(dict), map(GridSampling, elements(grid(S))))
    tensorproduct(duals...)
end


## Complexified dictionaries

import BasisFunctions: ComplexifiedDict

dualdictionary(dict::ComplexifiedDict) =
    ComplexifiedDict(dualdictionary(superdict(dict)))

dualsamplingoperator(dict::ComplexifiedDict, S) =
    dualsamplingoperator(superdict(dict), S)


## Extension frames

dualdictionary(dict::ExtensionFrame) = ExtensionFrame(support(dict),dualdictionary(basis(dict)))

function dualsamplingoperator(d::ExtensionFrame, S)
    grid1 = grid(dest(S))
    grid2 = supergrid(grid1)
    T = coefficienttype(dest(S))
    op = quadraturenormalization(T, grid2)
    if op isa ScalingOperator
        scalar(op) * S
    else
        E = extension_operator(GridBasis{T}(grid1),GridBasis{T}(grid2))
        R = E'
        tot = R*quadraturenormalization(T, grid2)*E
        tot * S
    end
end
