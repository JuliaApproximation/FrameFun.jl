# fourierdomains.jl

# Warning: currently not functional. TODO: update


"""
A FourierDomain is a domain that is implicitly defined by a Fourier series.
Example: if relop is x -> x-2 > 0, then the Fourier domain is the domain where the given
Fourier series is greater than 2.
"""
immutable FourierDomain{N,T} <: AbstractDomain{N}
    f       ::  ExpFun{N,T}
    relop	::	Function
    box     ::  FBox{N,T}
end

FourierDomain{N,T}(f::ExpFun{N,T}, relop::Function) = FourierDomain{N,T}(f, relop, dbox(f))

function in(x::AbstractVector, d::FourierDomain)
	z = call(d.f, x)
	d.relop(z)
end

function in(g::Grid, d::FourierDomain)
	z = call(d.f, g)[1]
	map(d.relop, z)
end

(<){T <: Number}(f::ExpFun, a::T) = FourierDomain(f, x-> x < a)
(<=){T <: Number}(f::ExpFun, a::T) = FourierDomain(f, x-> x <= a)
(>){T <: Number}(f::ExpFun, a::T) = FourierDomain(f, x-> x > a)
(>=){T <: Number}(f::ExpFun, a::T) = FourierDomain(f, x-> x >= a)



"""
A ComparisonDomain is a domain that is implicitly defined by two Fourier series.
Example: if binop is (x,y) -> x > y, then the Fourier domain is the domain where the first given
Fourier series is greater than the second.
"""
immutable ComparisonDomain{N,T} <: AbstractDomain{N}
    f       ::  ExpFun{N,T}
    g		::	ExpFun{N,T}
    binop	::	Function
    box     ::  FBox{N,T}
end

# TODO: We don't do any bound checking here. We should. Are the domains compatible?
ComparisonDomain{N,T}(f::ExpFun{N,T}, g::ExpFun{N,T}, binop::Function) = FourierDomain{N,T}(f, g, binop, dbox(f))

function in(x::AbstractVector, d::ComparisonDomain)
	z1 = call(d.f, x)
	z2 = call(d.g, x)
	d.binop(z1, z2)
end

function in(g::Grid, d::ComparisonDomain)
	z1 = call(d.f, g)
	z2 = call(d.g, g)
	map(d.binop, z1, z2)
end

(<)(f::ExpFun, g::ExpFun) = ComparisonDomain(f, g, (x,y)-> x < y)
(<=)(f::ExpFun, g::ExpFun) = ComparisonDomain(f, g, (x,y)-> x <= y)
(>)(f::ExpFun, g::ExpFun) = ComparisonDomain(f, g, (x,y)-> x > y)
(>=)(f::ExpFun, g::ExpFun) = ComparisonDomain(f, g, (x,y)-> x >= y)


