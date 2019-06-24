module DictFuns
using StaticArrays, BasisFunctions, DomainSets

import BasisFunctions: expansion, element, restrict
import Base: eltype, length, size, adjoint, *, show

export DictFun
"""
A `DictFun` corresponds to an expansion in a dictionary, but it adds a simple user
interface for computing with functions.
"""
struct DictFun{S,T} <: Function
    expansion   ::  Expansion
end

DictFun(e::Expansion, args...) = DictFun{domaintype(dictionary(e)),codomaintype(dictionary(e))}(e, args...)

DictFun(dict::Dictionary{S,T}, coefficients = zeros(dict), args...) where {S,T} =
    DictFun{S,T}(Expansion(dict, coefficients), args...)

# Warning: not all 2d function sets have SVector{2,T} type, they could have (S,T) type
const DictFun1d{S <: Number,T} = DictFun{S,T}
const DictFun2d{S <: Number,T} = DictFun{SVector{2,S},T}
const DictFun3d{S <: Number,T} = DictFun{SVector{3,S},T}
const DictFun4d{S <: Number,T} = DictFun{SVector{4,S},T}

expansion(fun::DictFun) = fun.expansion

Base.broadcast(fun::DictFun, x...) = broadcast(expansion(fun), x...)

element(fun::DictFun, i) = DictFun(element(fun.expansion, i))
superfun(fun::DictFun) = DictFun(superdict(dictionary(fun)), coefficients(fun))

for op in (:dictionary, :dimension, :coefficients, :eltype, :numtype, :length, :size, :support)
    @eval BasisFunctions.$op(fun::DictFun) = $op(fun.expansion)
end

for op in (:adjoint, :∫, :∂x, :∂y, :∂z, :∫∂x, :∫∂y, :∫∂z, :differentiate, :antidifferentiate)
    @eval BasisFunctions.$op(fun::DictFun{S,T}, args...) where {S,T} = DictFun{S,T}($op(fun.expansion, args...))
end

for op in (:+, :-, :*)
    @eval Base.$op(fun1::DictFun,fun2::DictFun) = DictFun($op(fun1.expansion,fun2.expansion))
end

for op in (:+, :-, :*)
    @eval Base.$op(a::Number,fun::DictFun) = DictFun($op(a,fun.expansion))
end

for op in (:+, :-, :*)
    @eval Base.$op(fun::DictFun,a::Number) = $op(a,fun)
end


# Delegate operator applications to the underlying expansion
function (*)(op::DictionaryOperator, fun::DictFun)
    @assert (src(op) == dictionary(fun)) | (src(op) == basis(dictionary(fun)))
    DictFun(dest(op),op*coefficients(fun))
end

# Delegate all calling to the underlying expansion.
(fun::DictFun)(x...; options...) = expansion(fun)(x...; options...)


# show(io::IO, fun::DictFun) = show(io, fun, set(fun))

function show(io::IO, fun::DictFun, set::Dictionary)
    println(io, "A ", dimension(fun), "-dimensional DictFun with ", length(coefficients(fun)), " degrees of freedom.")
    println(io, "Basis: ", name(set))
end

end
