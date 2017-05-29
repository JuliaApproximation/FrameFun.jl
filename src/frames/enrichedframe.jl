# enrichedframe.jl

"""
An EnrichedFrame is the concatenation of a basis with a bounded number of extra
functions. It differs from a SumFrame only in the way it resizes: only the basis
is resized, not the other set(s).
"""
struct EnrichedFrame{N,T} <: DerivedSet{N,T}
    superset    ::  FunctionSet{N,T}

    function EnrichedFrame{N,T}(set) where {N,T}
        @assert is_enrichedframe(set)
        new(set)
    end
end

EnrichedFrame{N,T}(set::FunctionSet{N,T}) = EnrichedFrame{N,T}(set)

is_enrichedframe(set::FunctionSet) = false
is_enrichedframe(set::MultiSet) = has_transform(element(set,1))


similar_set(f::EnrichedFrame, set::FunctionSet) = EnrichedFrame(set)

# The basis of the enriched frame is the first element
basis(f::EnrichedFrame) = element(superset(f),1)

name(f::EnrichedFrame) = "The basis " * name(element(f,1)) * " enriched with " * name(element(f,2))


############
# Extension
############

# We only extend the basis, not the other set(s)
extension_size(f::EnrichedFrame) = extension_size(basis(f))

resize(f::EnrichedFrame, n) =
    similar_set(f, MultiSet([resize(basis(f), n); elements(f)[2:end]]))

################
# Approximation
################

approximation_operator(frame::EnrichedFrame; options...) = SumFrameSolver(frame)
