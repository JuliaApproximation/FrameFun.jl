# fourierdomains.jl


"""
A FourierDomain is a domain that is implicitly defined by a Fourier series.
Example: if relop is x -> x-2 > 0, then the Fourier domain is the domain where the given
Fourier series is greater than 2.
"""
struct FourierDomain{T} <: Domain{T}
    f       ::  SetFun
    relop	::	Function
end

function indomain(x, d::FourierDomain)
	z = d.f(x)
	in(x,domain(d.f)) && d.relop(z)
end

# TODO: add postprocessing so that x outside domain -> NaN
function in(g::AbstractGrid, d::FourierDomain)
	z = d.f(g)[1]
	map(d.relop, z)
end

(<)(f::SetFun{T}, a::Number) where {T} = FourierDomain{T}(f, x-> x < a)
(<=)(f::SetFun{T}, a::Number) where {T} = FourierDomain{T}(f, x-> x <= a)
(>)(f::SetFun{T}, a::Number) where {T} = FourierDomain{T}(f, x-> x > a)
(>=)(f::SetFun{T}, a::Number) where {T} = FourierDomain{T}(f, x-> x >= a)

boundingbox(d::FourierDomain) = boundingbox(domain(d.f))

"""
A ComparisonDomain is a domain that is implicitly defined by two Fourier series.
Example: if binop is (x,y) -> x > y, then the Fourier domain is the domain where the first given
Fourier series is greater than the second.
"""
struct ComparisonDomain{T} <: Domain{T}
    f       ::  SetFun
    g		::	SetFun
    binop	::	Function
end

function indomain(x, d::ComparisonDomain)
	z1 = d.f(x)
	z2 = d.g(x)
	in(x,domain(d.f)) && in(x,domain(d.g)) && d.binop(z1, z2)
end

# TODO: add postprocessing so that x outside domain -> NaN
function in(g::AbstractGrid, d::ComparisonDomain)
	z1 = d.f(g)
	z2 = d.g(g)
	in(g,domain(d.f)) & in(g,domain(d.g)) & map(d.binop, z1, z2)
end

(<)(f::SetFun{T}, g::SetFun{T}) where {T} = ComparisonDomain{T}(f, g, (x,y)-> x < y)
(<=)(f::SetFun{T}, g::SetFun{T}) where {T} = ComparisonDomain{T}(f, g, (x,y)-> x <= y)
(>)(f::SetFun{T}, g::SetFun{T}) where {T} = ComparisonDomain{T}(f, g, (x,y)-> x > y)
(>=)(f::SetFun{T}, g::SetFun{T}) where {T} = ComparisonDomain{T}(f, g, (x,y)-> x >= y)

boundingbox(d::ComparisonDomain) = boundingbox(domain(d.f))+boundingbox(domain(d.g))
