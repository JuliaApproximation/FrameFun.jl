# weightedsumframe.jl

"""
A WeightedSumFrame is the concatenation of multiple sets, possibly weighted
except for the first one.
"""
immutable WeightedSumFrame{N,T} <: DerivedSet{N,T}
    superset    ::  FunctionSet{N,T}

    function WeightedSumFrame(set)
        @assert is_weightedsumset(set)
        new(set)
    end
end

WeightedSumFrame{N,T}(set::FunctionSet{N,T}) = WeightedSumFrame{N,T}(set)

is_weightedsumset(set::FunctionSet) = false

is_weightedsumset(set::MultiSet) = _is_weightedsumset(set, elements(set)...)

_is_weightedsumset(set::MultiSet, set::FunctionSet...) = false
_is_weightedsumset(set::MultiSet, set::WeightedSet...) = false
_is_weightedsumset(set::MultiSet, set1::FunctionSet, set2::WeightedSet) = true
_is_weightedsumset(set::MultiSet, set1::FunctionSet, set2::WeightedSet, set3::WeightedSet) = true

similar_set(f::WeightedSumFrame, set::FunctionSet) = WeightedSumFrame(set)

sumframe(sets::FunctionSet...) = WeightedSumFrame(multiset(sets...))
