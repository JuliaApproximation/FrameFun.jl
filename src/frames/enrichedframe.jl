# enrichedframe.jl

"""
An EnrichedFrame is the concatenation of a basis with a bounded number of extra
functions.
"""
immutable EnrichedFrame{N,T} <: DerivedSet{N,T}
    superset    ::  FunctionSet{N,T}
end

EnrichedFrame{N,T}(set::FunctionSet{N,T}) = EnrichedFrame{N,T}(set)

is_enrichedframe(set::FunctionSet) = false

is_enrichedframe(set::MultiSet) = _is_enrichedframe(set, elements(set)...)
is_enrichedframe(set::MultiSet, set1::FunctionSet, sets::FunctionSet...) =
    is_basis(set1) && has_transform(set1)


similar_set(f::EnrichedFrame, set::FunctionSet) = EnrichedFrame(set)

# The basis of the enriched frame is the first element
basis(f::EnrichedFrame) = element(superset(f),1)

############
# Extension
############

# We only extend the basis, not the other set(s)
extension_size(f::EnrichedFrame) = extension_size(basis(f))

resize(f::EnrichedFrame, n) =
    similar_set(f, replace(superset(f), 1, resize(basis(f),n)))
