# arithmetic.jl
# Functions related to computations with ExpFun's and EFun's

(+){T <: Number}(fl::EFun, a::T) = addfun(fl, a)
(+){T <: Number}(a::T, fl::EFun) = fl+a
(+)(a::EFun, b::EFun) = addfun(a, b)

(+){T <: Number}(f::ExpFun, a::T) = addfun(f, a)
(+){T <: Number}(a::T, f::ExpFun) = f+a
(+)(a::ExpFun, b::ExpFun) = addfun(a, b)

(*){T <: Number}(fl::EFun, a::T) = timesfun(fl, a)
(*){T <: Number}(a::T, fl::EFun) = fl*a
(*){T <: Number}(f::ExpFun, a::T) = timesfun(f, a)
(*){T <: Number}(a::T, f::ExpFun) = f*a


function addfun{T <: Number}(fl::EFun, a::T)
	layer = copy(fl)
	layer.coef[1] = layer.coef[1]+a
	layer
end

# We really should be resampling here in case the dimensions don't match!
# Right now this will just throw an error
addfun(a::EFun, b::EFun) = EFun(domain(a) & domain(b), a.param, copy(a.context), a.coef+b.coef)

timesfun{T <: Number}(fl::EFun, a::T) = similar(fl, fl.coef*a)

addfun{N,T,S <: Number}(f::ExpFun{N,T}, a::S) = maplayers(f, layer -> layer+a)

function addfun{N,T}(a::ExpFun{N,T}, b::ExpFun{N,T})
	g = emptyfun(N, T)
	for i = 1:nblayers(a)
		push!(g, layer(a,i)+layer(b,i))
	end
	g
end

timesfun{N,T,S <: Number}(f::ExpFun{N,T}, a::S) = maplayers(f, layer -> layer*a)

differentiate(fl::EFun, i = 1) = similar(fl, differentiate(fl.coef, i))

differentiate{N,T}(f::ExpFun{N,T}, i = 1) = maplayers(fl, layer -> differentiate(layer, i))


setcoef!(fl::EFun, coef) = fl.coef[:] = coef



getindex{N,T}(f::ExpFun{N,T}, d::AbstractDomain) = maplayers(f, layer -> similar(layer, layer.domain & d))

function setindex!(f::ExpFun, g::Function, d::AbstractDomain)
	layer = computelayer(g, d)
	push!(f, layer)
	g
end


function setindex!(f::ExpFun, g::ExpFun, d::AbstractDomain)
	for i = 1:nblayers(g)
		layer = g.layers[i]
		layer2 = similar(layer, layer.domain & d)
		push!(f, layer2)
	end
	g
end


function setindex!{T <: Number}(f::ExpFun, a::T, d::AbstractDomain)
	dof, extension, oversampling, method = default_parameters()
	fp,fc = paramcontext(d, dof, extension, oversampling)
	layer = EFun(d, fp, fc, zeros(Complex{T}, fp.n...))
	layer.coef[1] = a
	push!(f, layer)
	a
end


