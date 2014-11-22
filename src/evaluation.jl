# evaluation.jl
# Routines to evaluate ExpFun's

function relocate_point{N,T}(x, b::FBox{N,T})
	T[(x[i]-left(b,i)) / size(b,i) for i=1:N]
end

##$ Calling a function

## Without arguments: return nothing

call(f::AbstractFun) = nothing

## With a single value as argument: evaluate in the point x
call{T <: Number}(f::AbstractFun{1,T}, x::T) = call(f, [x,])

function call(fl::EFun, x::AbstractVector)
	box = dbox(fl)
	u = relocate_point(x, box)
	z = evaluate(fl.param, fl.context, fl.coef, u)
end


function call(f::ExpFun, x::AbstractVector)
	for i = nblayers(f):-1:1
		fl = layer(f, i)
		if in(x, fl.domain)
			z = call(fl, x)
			return real(z)
		end
	end
	NaN*one(T)
end


## With a grid or a tuple as argument: evaluate on the (corresponding) grid

# Compute the number of equispaced points on [a,c] that includes n points on [a,b].
numberofpoints(a, b, c, n) = (h = (b-a)/(n-1); int(ceil( (c-a) / h )))

# Compute suitable evaluation grids with m points inside the dbox
function suitablegrids(dbox::FBox, ebox::FBox, m)
	l = tuple([ numberofpoints(left(dbox,i), right(dbox,i), right(ebox,i), m[i]) for i=1:length(m) ]...)
	gridl = periodicgrid(ebox, l)
	gridm = subgrid(gridl, m)
	(gridm,gridl)
	# exact_match = reduce( (&&), [right(gridn,j)==right(dbox,j) for j=1:dim(dbox)])
end


suitablegrids(f::AbstractFun, m::NTuple) = suitablegrids(dbox(f), ebox(f), m)
# suitablegrids(fl::EFun, m::NTuple) = suitablegrids(dbox(fl), ebox(fl), m)
# suitablegrids(f::ExpFun, m::NTuple) = suitablegrids(dbox(f), ebox(f), m)

suitablegrids(fl::ExpFun) = suitablegrids(fl, fl.param.m)
suitablegrids(f::ExpFun) = suitablegrids(f.layer[1])

# Special case for 1D
suitablegrids(f::AbstractFun{1}, m::Int) = suitablegrids(f, (m,))


# Evaluate the Fourier extension on the extended domain
eval_extension(fl::EFun, l) = evaluate(fl.param, fl.context, fl.coef, l)

# Evaluate just the first layer when asked
eval_extension(f::ExpFun, l::NTuple) = eval_extension(layer(f,1), l)

# Return a mask for the given grid, use cached one if available
mask(fl::EFun, grid::Grid) = (grid == Grid(dbox(fl), fl.param.m) ? fl.context.maskM : in(grid, fl.domain))


# Evaluate the Fourier series in a given number of points on the bounding box.
# A suitable grid is computed and returned as second output.
function evaluate_unmasked{N,T}(fl::EFun{N,T}, m::NTuple)
	gridm,gridl = suitablegrids(fl, m)
	l = size(gridl)
	z1 = eval_extension(fl, l)
	z = zeros(Complex{T}, m...)
	truncate!(z, z1)
	z,gridm
end

# The function call evaluates the Fourier series in the given grid. FFT is used if possible,
# point by point evaluation otherwise. Points not in the domain evaluate to NaN.
function evaluate_unmasked(fl::EFun, grid::Grid)
	m = size(grid)
	gridm,gridl = suitablegrids(fl, m)
	if grid == gridm
		# We can use the fft routine
		z = evaluate_unmasked(fl, m)[1]
	else
		# The grid is not suitable for fft, we evaluate point by point
		z = evaluate_in_points(fl, grid)
	end
	z
end

function evaluate_masked(fl::EFun, m::NTuple)
	z,grid = evaluate_unmasked(fl, m)
	mask = mask(fl, grid)
	z[~mask] = NaN
end

function evaluate_masked(fl::EFun, grid::Grid)
	z = evaluate_unmasked(fl, grid)
	mask = mask(fl, grid)
	z[~mask] = NaN
end


call(fl::EFun, g::Grid) = evaluate_masked(fl, g)


@ngenerate N Array{Complex{T},N} function evaluate_in_points{N,T}(fl::EFun{N,T}, g::Grid{N,T})
	l = size(g)
	z = zeros(Complex{T}, l...)
	@nloops N i g begin
		(@nref N z i) = call(fl, (@nref N g i))
	end
	z
end

call(f::ExpFun, g::Grid) = evaluate(f, g)

function evaluate{N,T}(f::ExpFun{N,T}, grid::Grid)
	m = size(grid)

	z = fill(convert(Complex{T},NaN), m...)
	for i = 1:nblayers(f)
		z1 = evaluate_unmasked(layer(f,i), grid)
		maskm = mask(layer(f, i), grid)
		z[maskm] = z1[maskm]
	end
	real(z)
end

function evaluate{N,T}(f::ExpFun{N,T}, m::NTuple{N})
	z = fill(convert(Complex{T},NaN), m...)
	gridm = Grid(dbox(f), m)	# just define this variable for now, if it is only defined in the loop it is not known outside (and can't be returned)
	for i = 1:nblayers(f)
		z1,grid = evaluate_unmasked(layer(f,i), m)
		maskm = mask(layer(f, i), grid)
		z[maskm] = z1[maskm]
		gridm = grid
	end
	real(z),gridm
end
# evaluate(f::ExpFun, m::NTuple) = call(f, Grid(dbox(f), m))
