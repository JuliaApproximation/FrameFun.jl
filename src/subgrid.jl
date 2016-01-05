# subgrid.jl


abstract AbstractSubGrid{N,T,G} <: AbstractGrid{N,T}

eltype{N,T,G}(::Type{AbstractSubGrid{N,T,G}}) = eltype(G)

grid(g::AbstractSubGrid) = g.grid


"""
A MaskedGrid is a subgrid of another grid that is defined by a mask.
The mask is true or false for each point in the supergrid. The set of points
for which it is true make up the MaskedGrid.
"""
immutable MaskedGrid{G,ID,N,T} <: AbstractSubGrid{N,T,G}
    grid	::	G
    mask	::	Array{Bool,ID}
    indices ::  Vector{Vec{N,Int}}
    M		::	Int				# Total number of points in the mask

    MaskedGrid(grid::AbstractGrid{N,T},  mask, indices) = new(grid, mask, indices, sum(mask))
end

function MaskedGrid{N,T}(grid::AbstractGrid{N,T}, mask, indices)
	@assert size(grid) == size(mask)

	MaskedGrid{typeof(grid),index_dim(grid),N,T}(grid, mask, indices)
end

# TODO: make this more elegant and general
# These are for the assignment to indices in the function below.
convert(::Type{Tuple{Int}}, i::CartesianIndex{1}) = (i[1],)
convert(::Type{Tuple{Int,Int}}, i::CartesianIndex{2}) = (i[1],i[2])
convert(::Type{Tuple{Int,Int,Int}}, i::CartesianIndex{3}) = (i[1],i[2],i[3])
convert(::Type{Tuple{Int,Int,Int,Int}}, i::CartesianIndex{4}) = (i[1],i[2],i[3],i[4])

function MaskedGrid{N,T}(grid::AbstractGrid{N,T}, domain::AbstractDomain{N,T})
    mask = in(grid, domain)
    indices = Array(Vec{N,Int}, sum(mask))
    i = 1
    for m in eachindex(grid)
        if mask[m]
            indices[i] = m
            i += 1
        end
    end
    MaskedGrid(grid, mask, indices)
end


length(g::MaskedGrid) = g.M

size(g::MaskedGrid) = (length(g),)


# Check whether element grid[i] (of the underlying grid) is in the masked grid.
in(i, g::MaskedGrid) = g.mask[i]

getindex(g::MaskedGrid, idx::Int) = getindex(g.grid, g.indices[idx]...)


# Efficient extension operator
function apply!{G <: MaskedGrid}(op::Extension, dest, src::DiscreteGridSpace{G}, coef_dest::AbstractArray, coef_src::AbstractArray)
    @assert length(coef_src) == length(src)
    @assert length(coef_dest) == length(dest)
    @assert grid(dest) == grid(grid(src))

    grid1 = grid(src)
    fill!(coef_dest, 0)

    l = 0
    for i in eachindex(grid1.grid)
        if in(i, grid1)
            l += 1
            coef_dest[i] = coef_src[l]
        end
    end
end


# Efficient restriction operator
function apply!{G <: MaskedGrid}(op::Restriction, dest::DiscreteGridSpace{G}, src, coef_dest::AbstractArray, coef_src::AbstractArray)
    @assert length(coef_src) == length(src)
    @assert length(coef_dest) == length(dest)
    @assert grid(src) == grid(grid(dest))

    grid1 = grid(dest)

    l = 0
    for i in eachindex(grid1.grid)
        if in(i, grid1)
            l += 1
            coef_dest[l] = coef_src[i]
        end
    end
end



"""
An IndexSubGrid is a subgrid corresponding to a certain range of indices of the
underlying (one-dimensional) grid.
"""
immutable IndexSubGrid{G,T} <: AbstractSubGrid{1,T,G}
	grid	::	G
	i1		::	Int
	i2		::	Int

	function IndexSubGrid(grid::AbstractGrid1d{T}, i1, i2)
		@assert 1 <= i1 <= length(grid)
		@assert 1 <= i2 <= length(grid)
		@assert i1 <= i2

		new(grid, i1, i2)
	end
end

IndexSubGrid{T}(grid::AbstractGrid1d{T}, i1, i2) = IndexSubGrid{typeof(grid), T}(grid, i1, i2)

left(g::IndexSubGrid) = g.grid[g.i1]

right(g::IndexSubGrid) = g.grid[g.i2]

length(g::IndexSubGrid) = g.i2 - g.i1 + 1

getindex(g::IndexSubGrid, i) = g.grid[g.i1+i-1]

# Check whether element grid[i] (of the underlying grid) is in the indexed subgrid.
in(i, g::IndexSubGrid) = (i >= g.i1) && (i <= g.i2)

stepsize{G <: AbstractEquispacedGrid}(g::IndexSubGrid{G}) = stepsize(g.grid)

range{G <: AbstractEquispacedGrid}(g::IndexSubGrid{G}) = left(g) : stepsize(g) : right(g)


# Efficient extension operator
function apply!{G <: IndexSubGrid}(op::Extension, dest::DiscreteGridSpace, src::DiscreteGridSpace{G}, coef_dest::AbstractArray, coef_src::AbstractArray)
    @assert length(coef_src) == length(src)
    @assert length(coef_dest) == length(dest)
    @assert grid(dest) == grid(grid(src))

    grid1 = grid(src)
    fill!(coef_dest, 0)

    l = 0
    for i in grid1.i1:grid1.i2
        l += 1
        coef_dest[i] = coef_src[l]
    end
end


# Efficient restriction operator
function apply!{G <: IndexSubGrid}(op::Restriction, dest::DiscreteGridSpace{G}, src::DiscreteGridSpace, coef_dest::AbstractArray, coef_src::AbstractArray)
    @assert length(coef_src) == length(src)
    @assert length(coef_dest) == length(dest)
    @assert grid(src) == grid(grid(dest))

    grid1 = grid(dest)

    l = 0
    for i in grid1.i1:grid1.i2
        l += 1
        coef_dest[l] = coef_src[i]
    end
end



"Create a suitable subgrid that covers a given domain."
function subgrid(grid::AbstractEquispacedGrid, domain::Interval)
    a = left(domain)
    b = right(domain)
    h = stepsize(grid)
    idx_a = convert(Int, ceil( (a-left(grid))/stepsize(grid))+1 )
    idx_b = convert(Int, floor( (b-left(grid))/stepsize(grid))+1 )
    IndexSubGrid(grid, idx_a, idx_b)
end

subgrid(grid::AbstractGrid, domain::AbstractDomain) = MaskedGrid(grid, domain)


