
"Supertype of expansions in frames and bases."
abstract type AbstractFun end

dictionary(fun::AbstractFun) = dictionary(expansion(fun))

(fun::AbstractFun)(x...) = expansion(fun)(x...)

"The combination of an expansion with a platform."
struct SimpleFun <: AbstractFun
    platform
    expansion
end

platform(fun::SimpleFun) = fun.platform
expansion(fun::SimpleFun) = fun.expansion

function simplefun(f, platform, args...; options...)
    F = Fun(f, platform, args...; options...)
    SimpleFun(platform, F)
end

similarfun(f::SimpleFun, F::Expansion) = SimpleFun(platform(f), F)

funvariable(platform::Platform) = simplefun(x->x, platform)

chebvariable(::Type{T} = Float64) where {T} =
    funvariable(platform(ChebyshevT{T}(2)))

## Arithmetics

for op in (:+, :-, :*)
    # combination with number
    @eval Base.$op(a::Number, fun::AbstractFun) = similarfun(fun, $op(a, expansion(fun)))
    @eval Base.$op(fun::AbstractFun, a::Number) = similarfun(fun, $op(expansion(fun), a))
    # combination of funs
    @eval Base.$op(f::AbstractFun, g::AbstractFun) = similarfun(f, $op(expansion(f), expansion(g)))
end
Base.:/(fun::AbstractFun, a::Number) = similarfun(fun, expansion(fun)/a)

import Base: cos, sin, tan
for op in [:cos, :sin, :tan]
    @eval function $op(f::AbstractFun)
        F = Fun(x->cos(f(x)), platform(f))
        similarfun(f, F)
    end
end
