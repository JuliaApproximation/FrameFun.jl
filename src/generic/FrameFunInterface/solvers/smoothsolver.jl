# Remove this part?
export sobolevAZweight
sobolevAZweight(A::DictionaryOperator; options...) = inv(WeightedSmoothingOperator(src(A); options...))
# For Dictionary's that have a DC component
dc_index(b::ChebyshevT) = 1
dc_index(b::Fourier) = 1
function WeightedSmoothingOperator(dict::Dictionary; scaling = default_scaling_function, order = 1, options...)
    if order == 1
        WeightedSmoothingOperator(dict, scaling)
    else
        WeightedSmoothingOperator(dict, (d,i) -> scaling(d,i)^order)
    end
end
WeightedSmoothingOperator(dict::Dictionary, scaling) = DiagonalOperator(dict, dict, [scaling(dict, i) for i in eachindex(dict)][:])

function default_scaling_function(dict::Dictionary1d, idx)
    f = abs(value(native_index(dict, idx)))
    1 + f
end
default_scaling_function(dict::TensorProductDict, I) =
    default_scaling_function(element(dict, 1), I[1]) + default_scaling_function(element(dict, 2), I[2])
