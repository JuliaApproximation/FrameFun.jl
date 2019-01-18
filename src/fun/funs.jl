
abstract type AbstractFun
end

"""
A `DictFun` corresponds to an expansion in a dictionary, but it adds a simple user
interface for computing with functions.
"""
struct DictFun{S,T} <: AbstractFun
    expansion   ::  Expansion
end

DictFun(e::Expansion, args...) = DictFun{domaintype(dictionary(e)),codomaintype(dictionary(e))}(e, args...)

DictFun(dict::Dictionary{S,T}, coefficients = zeros(dict), args...) where {S,T} =
    DictFun{S,T}(Expansion(dict, coefficients), args...)

DictFun(domain::Domain, basis::Dictionary, args...) = DictFun(ExtensionFrame(domain, basis), args...)

# Warning: not all 2d function sets have SVector{2,T} type, they could have (S,T) type
const DictFun1d{S <: Number,T} = DictFun{S,T}
const DictFun2d{S <: Number,T} = DictFun{SVector{2,S},T}
const DictFun3d{S <: Number,T} = DictFun{SVector{3,S},T}
const DictFun4d{S <: Number,T} = DictFun{SVector{4,S},T}

expansion(fun::DictFun) = fun.expansion

Base.broadcast(fun::DictFun, x...) = broadcast(expansion(fun), x...)

element(fun::DictFun, i) = DictFun(element(fun.expansion, i))
superfun(fun::DictFun) = DictFun(superdict(dictionary(fun)), coefficients(fun))

for op in (:dictionary, :dimension, :coefficients, :eltype, :numtype, :length, :size)
    @eval $op(fun::DictFun) = $op(fun.expansion)
end

for op in (:adjoint, :∫, :∂x, :∂y, :∂z, :∫∂x, :∫∂y, :∫∂z, :differentiate, :antidifferentiate)
    @eval $op(fun::DictFun{S,T}, args...) where {S,T} = DictFun{S,T}($op(fun.expansion, args...))
end

for op in (:ExtensionFrame, :basis, :basisspan)
    @eval $op(fun::DictFun) = $op(fun, dictionary(fun))
end

domain(fun::DictFun) = domain(fun, dictionary(fun))

for op in (:+, :-, :*)
    @eval $op(fun1::DictFun,fun2::DictFun) = DictFun($op(fun1.expansion,fun2.expansion))
end

for op in (:+, :-, :*)
    @eval $op(a::Number,fun::DictFun) = DictFun($op(a,fun.expansion))
end

for op in (:+, :-, :*)
    @eval $op(fun::DictFun,a::Number) = $op(a,fun)
end

ExtensionFrame(fun::DictFun, set::ExtensionFrame) = set

basis(fun::DictFun, set::ExtensionFrame) = basis(set)

basis(fun::DictFun, dict::Dictionary) = dict

domain(fun::DictFun, dict::Dictionary) = support(dict)

# Delegate operator applications to the underlying expansion
function (*)(op::DictionaryOperator, fun::DictFun)
    @assert (src(op) == dictionary(fun)) | (src(op) == basis(dictionary(fun)))
    DictFun(dest(op),op*coefficients(fun))
end

# Delegate all calling to the underlying expansion.
(fun::DictFun)(x...) = expansion(fun)(x...)


# show(io::IO, fun::DictFun) = show(io, fun, set(fun))

function show(io::IO, fun::DictFun, set::Dictionary)
    println(io, "A ", dimension(fun), "-dimensional DictFun with ", length(coefficients(fun)), " degrees of freedom.")
    println(io, "Basis: ", name(set))
end

function show(io::IO, fun::DictFun, set::ExtensionFrame)
    println(io, "A ", dimension(fun), "-dimensional DictFun with ", length(coefficients(fun)), " degrees of freedom.")
    println(io, "Basis: ", name(basis(set)))
    println(io, "Domain: ", support(set))
end

getindex(expansion::Expansion, domain::Domain) = restrict(expansion, domain)

getindex(fun::DictFun, domain::Domain) = restrict(expansion(fun), domain)

restrict(expansion::Expansion, domain::Domain) = _restrict(expansion, dictionary(expansion), domain)

function _restrict(expansion::Expansion, set::ExtensionFrame, domain1::Domain)
    @assert dimension(set) == dimension(domain1)

    domain2 = support(set)
    newdomain = domain1 ∩ domain2
    DictFun(newdomain, basis(set), coefficients(expansion))
end

function _restrict(expansion::Expansion, set::Dictionary, domain::Domain)
    @assert dimension(set) == dimension(domain)
    # We should check here whether the given domain lies in the support of the set
    DictFun(domain, set, coefficients(expansion))
end
