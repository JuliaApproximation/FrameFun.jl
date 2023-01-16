
"Supertype of expansions in frames and bases."
abstract type AbstractFun end

dictionary(fun::AbstractFun) = dictionary(expansion(fun))

(fun::AbstractFun)(x::Number) = expansion(fun)(x)
(fun::AbstractFun)(x::Vararg{Number,N}) where {N} = expansion(fun)(SVector{N}(x))
(fun::AbstractFun)(x::AbstractVector) = expansion(fun)(x)

"The combination of an expansion with a platform."
struct SimpleFun <: AbstractFun
    platform
    expansion
end

platform(fun::SimpleFun) = fun.platform
expansion(fun::SimpleFun) = fun.expansion

SimpleFun(e::Expansion) = SimpleFun(platform(dictionary(e)), e)

function simplefun(f, platform, args...; options...)
    F = Fun(f, platform, args...; options...)
    SimpleFun(platform, F)
end

function similarfun(f::SimpleFun, F::Expansion)
    F2 = expansion(dictionary(F), coefficients(F))
    if !isreal(F2)
        SimpleFun(complex(platform(f)), F2)
    else
        SimpleFun(platform(f), F2)
    end
end

funvariable(platform::Platform) = simplefun(x->x, platform)

function funvariables(platform::Platform, dim::Int)
    vars = simplefun.([ (x...) -> x[k] for k in 1:dim], Ref(platform))
    tuple(vars...)
end

chebvariable(::Type{T} = Float64) where {T} =
    SimpleFun(BasisFunctions.expansion_of_x(ChebyshevT{T}(2)))

function chebvariable(domain::AbstractInterval{T}) where {T}
    x = BasisFunctions.expansion_of_x(ChebyshevT{float(T)}(2) → domain)
    SimpleFun(x)
end

function chebvariables(dim::Int, ::Type{T} = Float64) where {T}
    basis_1 = ChebyshevT{T}(1)
    basis_x = ChebyshevT{T}(2)
    c_1 = BasisFunctions.coefficients_of_one(basis_1)
    c_x = BasisFunctions.coefficients_of_x(basis_x)
    if dim == 2
        c2_x = collect(OuterProductArray(c_x, c_1))
        x = Expansion(basis_x ⊗ basis_1, c2_x)
        c2_y = collect(OuterProductArray(c_1, c_x))
        y = Expansion(basis_1 ⊗ basis_x, c2_y)
        SimpleFun(x), SimpleFun(y)
    else
        error("Not implemented")
    end
end


## Arithmetics

import BasisFunctions: roots
roots(fun::AbstractFun) = roots(expansion(fun))

for op in (:adjoint, :complex)
    @eval Base.$op(fun::AbstractFun) = similarfun(fun, $op(expansion(fun)))
end

for op in (:+, :-, :*)
    # combination with number
    @eval Base.$op(a::Number, fun::AbstractFun) = similarfun(fun, $op(a, expansion(fun)))
    @eval Base.$op(fun::AbstractFun, a::Number) = similarfun(fun, $op(expansion(fun), a))
    # combination of funs
    @eval Base.$op(f::AbstractFun, g::AbstractFun) = similarfun(f, $op(expansion(f), expansion(g)))
end
Base.:/(fun::AbstractFun, a::Number) = similarfun(fun, expansion(fun)/a)

"For a given fun `f(x)` compute the fun `g(f(x))`."
function compute_fun(fun::AbstractFun, g)
    F = Fun((x...)->g(fun(x...)), platform(fun))
    similarfun(fun, F)
end

Base.:^(fun::AbstractFun, a::Integer) = compute_fun(fun, x->x^a)
Base.:/(a::Number, fun::AbstractFun) = compute_fun(fun, x->1/x)

# list of functions to override was inspired by ApproxFun
for op in [:sqrt, :cbrt, :abs2, :inv, :log, :log10, :log2, :log1p,
            :exp, :exp2, :expm1, :sin, :cos, :tan, :sec, :csc, :cot,
            :sind, :cosd, :tand, :secd, :cscd, :cotd, :asin, :acos, :atan,
            :asec, :acsc, :acot, :asind, :acosd, :atand, :asecd, :acscd,
            :acotd, :sinh, :cosh, :tanh, :sech, :csch, :coth, :asinh,
            :acosh, :atanh, :asech, :acsch, :acoth, :deg2rad, :rad2deg]
    @eval Base.$op(fun::AbstractFun) = compute_fun(fun, $op)
end

# list of functions to override was inspired by ApproxFun
for op in [:erf, :erfinv, :erfc, :erfcinv, :erfi, :gamma, :lgamma, :digamma,
        :invdigamma, :trigamma, :airyai, :airybi, :airyaiprime, :airybiprime,
        :besselj0, :besselj1, :bessely0, :bessely1, :erfcx, :dawson]
    @eval SpecialFunctions.$op(fun::AbstractFun) = compute_fun(fun, $op)
end
