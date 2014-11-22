module FrameFun

using ArrayViews
using ImmutableArrays
using Base.Cartesian
# using Devectorize

# These functions will be extended:
import Base.size
import Base.show
import Base.push!
import Base.getindex
import Base.eps
import Base.evaluate
import Base.in
import Base.copy
import Base.similar
import Base.string

export ExpFun, EFun, plot, plot_extension, plot_accuracy, plotslice, evaluate, call, suitablegrids

# Export domains
export Interval, Circle, Rectangle, Sphere, Cylinder, Mandelbrot, JuliaSet, Cube

export randomcircles, atomium, unitcircle, unitsquare, unitsphere, rotate

# Define eps for complex types as eps for the corresonding real type
eps{T}(::Type{Complex{T}}) = eps(T)

close{T}(a::T, b::T) = abs(a-b) < 10*eps(T)

include("my_lsqr.jl")

include("box.jl")

include("grid.jl")

include("domains.jl")

include("fractals.jl")

include("paramcontext.jl")

include("fourierextension.jl")



###############################################################################################
# Some common properties of Function objects are grouped in AbstractFun
###############################################################################################
abstract AbstractFun{N,T}

dim{N,T}(f::AbstractFun{N,T}) = N
numtype{N,T}(fl::AbstractFun{N,T}) = T

###############################################################################################
# An ExpFun consists of several EFun's
###############################################################################################

immutable EFun{N,T <: FloatingPoint,D} <: AbstractFun{N,T}
	domain  ::	D
	param   ::	FParam{N,T}
	context ::	FContext{N,T}
	coef    ::	Array{Complex{T},N}		# the Fourier coefficients
end


###############################################################################################
# An ExpFun has several layers. It has an overal bounding box and extension box, and all
# Fourier series in all of the layers are relative to these boxes.
###############################################################################################

type ExpFun{N,T <: FloatingPoint} <: AbstractFun{N,T}
	layers	::	Vector{EFun{N,T}}	# the different layers
	dbox	::	FBox{N,T}				# the overall bounding box
	ebox	::	FBox{N,T}				# the overall extension box
end

# Default parameters for the computation of new extensions
const DOF = 21
const EXTENSION = 2
const OVERSAMPLING = 2
const METHOD = "direct"


function ExpFun{N,T}(f, d::AbstractDomain{N,T}; dof = DOF, extension = EXTENSION*one(T), oversampling = OVERSAMPLING*one(T), method = METHOD)
	if N==3 && method == "direct" && dof == 21
		print("Warning: reducing dofs for this 3D problem.")
		dof = 5
	end
	layer = computelayer(f, d, dof, extension, oversampling, method)
	ExpFun{N,T}([layer], dbox(layer), ebox(layer))
end


emptyfun(N,T) = ExpFun{N,T}(Array(EFun{N,T},0), emptybox(N,T), emptybox(N,T))

nblayers(f::ExpFun) = length(f.layers)

layer(f::ExpFun, j) = f.layers[j]
layer(f::ExpFun) = f.layers

lastlayer(f::ExpFun) = layer(f, nblayers(f))

dbox(fl::EFun) = fl.context.dbox
ebox(fl::EFun) = fl.context.ebox
dbox(f::ExpFun) = f.dbox
ebox(f::ExpFun) = f.ebox

domain(fl::EFun) = fl.domain
domain(f::ExpFun, j) = domain(layer(f,j))
domain(f::ExpFun) = domain(f, nblayers(f))		# this should be generalized a little, the domain may be composite

ndof(fl::EFun) = prod(size(fl.coef))
ndof(f::ExpFun) = sum([ndof(layer(f,i)) for i=1:nblayers(f)])

size(fl::EFun) = size(fl.coef)
size(fl::EFun, j) = size(fl.coef, j)

copy(fl::EFun) = EFun(fl.domain, fl.param, copy(fl.context), copy(fl.coef))

# Create a new EFun with a new context and the given coefficients
similar(fl::EFun, coef::Array) = EFun(fl.domain, fl.param, copy(fl.context), coef)

# Create a new EFun with a new context and the given domain (but copy the coefficients)
similar(fl::EFun, d::AbstractDomain) = EFun(d, fl.param, similar(fl.context, d), copy(fl.coef))


# Create a new ExpFun by applying the function 'layermap' to each layer of f
function maplayers{N,T}(f::ExpFun{N,T}, layermap)
	g = emptyfun(N,T)
	for i = 1:nblayers(f)
		push!(g, layermap( layer(f,i)))
	end
	g
end

copy{N,T}(f::ExpFun{N,T}) = maplayers(f, layer -> copy(layer))


function computelayer(f, domain, dof, extension, oversampling, method)
	fp,fc = fe_paramcontext(domain, dof, extension, oversampling)
	if method == "direct"
		(c,A,B) = fe_approximate_direct(f, fp, fc)
	elseif method == "lsqr"
		c = fe_approximate_lsqr(f, fp, fc)
	elseif method == "extend"
		c = fe_approximate_extend(f, fp, fc)
	else
		print("Unknown method specified. Using iterative method instead.")
		c = fe_approximate_lsqr(f, fp, fc)
	end
	EFun(domain, fp, fc, c)
end


function recomputemasks!(fl::EFun)
	fp = fl.param
	fc = fl.context
	gridN = periodicgrid(fc.ebox, fp.n)
	gridM = Grid(fc.dbox, fp.m)
	gridL = periodicgrid(fc.ebox, fp.l)
	d = fc.d
	maskN = eval(fc.d, fc.gridN)
	maskM = eval(fc.d, fc.gridM)
	maskL = eval(fc.dd, fc.gridL)
	fc.maskN = maskN
	fc.maskM = maskM
	fc.maskL = maskL
end


function push!(f::ExpFun, layer::EFun)
	if nblayers(f) == 0
		# We must initialize the overall bounding box and extension box
		f.dbox = dbox(layer)
		f.ebox = ebox(layer)
	end
	assert(dbox(layer) == dbox(f))
	assert(ebox(layer) == ebox(f))
	push!(f.layers, layer)
end

function show{N}(io::IO, f::ExpFun{N})
	if nblayers(f) > 0
		print(io, "A ", N, "-dimensional expfun with ", nblayers(f), " layer(s) and a total of ", ndof(f), " degrees of freedom.\n")
		print(io, "The domain is: ", domain(f, 1), ".\n")
		print(io, "Thank you.")
	else
		print(io, "An empty ", N, "-dimensional expfun.")
	end
end

include("evaluation.jl")

include("arithmetic.jl")

include("fourierdomains.jl")

include("plotextensions.jl")

include("dpss.jl")

end # module FrameFun

