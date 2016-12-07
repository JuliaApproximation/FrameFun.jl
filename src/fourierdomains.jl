# fourierdomains.jl


"""
A FourierDomain is a domain that is implicitly defined by a Fourier series.
Example: if relop is x -> x-2 > 0, then the Fourier domain is the domain where the given
Fourier series is greater than 2.
"""
immutable FourierDomain{N,T} <: AbstractDomain{N}
    f       ::  FrameFun{N,T}
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

(<){T <: Number}(f::FrameFun, a::T) = FourierDomain(f, x-> x < a)
(<=){T <: Number}(f::FrameFun, a::T) = FourierDomain(f, x-> x <= a)
(>){T <: Number}(f::FrameFun, a::T) = FourierDomain(f, x-> x > a)
(>=){T <: Number}(f::FrameFun, a::T) = FourierDomain(f, x-> x >= a)

boundingbox(d::FourierDomain) = boundingbox(domain(d.f))

"""
A ComparisonDomain is a domain that is implicitly defined by two Fourier series.
Example: if binop is (x,y) -> x > y, then the Fourier domain is the domain where the first given
Fourier series is greater than the second.
"""
immutable ComparisonDomain{N,T} <: AbstractDomain{N}
    f       ::  FrameFun{N,T}
    g		::	FrameFun{N,T}
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

(<)(f::FrameFun, g::FrameFun) = ComparisonDomain(f, g, (x,y)-> x < y)
(<=)(f::FrameFun, g::FrameFun) = ComparisonDomain(f, g, (x,y)-> x <= y)
(>)(f::FrameFun, g::FrameFun) = ComparisonDomain(f, g, (x,y)-> x > y)
(>=)(f::FrameFun, g::FrameFun) = ComparisonDomain(f, g, (x,y)-> x >= y)

boundingbox(d::ComparisonDomain) = boundingbox(domain(d.f))+boundingbox(domain(d.g))
