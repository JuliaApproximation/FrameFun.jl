# subgrid.jl


abstract AbstractSubGrid{N,T} <: AbstractGrid{N,T}

dim{N,T}(::AbstractSubGrid{N,T}) = N
dim{N,T}(::Type{AbstractSubGrid{N,T}}) = N
dim{G <: AbstractSubGrid}(::Type{G}) = dim(super(G))

numtype{N,T}(::AbstractSubGrid{N,T}) = T
numtype{N,T}(::Type{AbstractSubGrid{N,T}}) = T
numtype{G <: AbstractSubGrid}(::Type{G}) = numtype(super(G))

# Default index dimension is 1
index_dim{N,T}(::AbstractSubGrid{N,T}) = 1
index_dim{N,T}(::Type{AbstractSubGrid{N,T}}) = 1
index_dim{G <: AbstractSubGrid}(::Type{G}) = 1

grid(g::AbstractSubGrid) = g.grid


"""
A MaskedGrid is a subgrid of another grid that is defined by a mask.
The mask is true or false for each point in the supergrid. The set of points
for which it is true make up the MaskedGrid.
"""
immutable MaskedGrid{G <: AbstractGrid,ID,N,T} <: AbstractSubGrid{N,T}
    grid	::	G
    mask	::	Array{Bool,ID}
    indices ::  Array{Int,2}
    M		::	Int				# Total number of points in the mask

    MaskedGrid(grid::AbstractGrid{N,T},  mask, indices) = new(grid, mask, indices, sum(mask))
end

function MaskedGrid{N,T}(grid::AbstractGrid{N,T}, mask, indices)
	@assert size(grid) == size(mask)

	MaskedGrid{typeof(grid),index_dim(grid),N,T}(grid, mask, indices)
end

function MaskedGrid{N,T}(grid::AbstractGrid{N,T}, domain::AbstractDomain{N,T})
    mask = in(grid, domain)
    indices = Array(Int,sum(mask),ndims(mask))
    i = 1
    for m in eachindex(mask)
        if mask[m]
            indices[i,:] = [ind2sub(mask,m)...]
            i += 1
        end
    end
    MaskedGrid(grid, mask, indices)
end


length(g::MaskedGrid) = g.M

size(g::MaskedGrid) = (length(g),)


# Check whether element grid[i] (of the underlying grid) is in the masked grid.
in(i, g::MaskedGrid) = g.mask[i]

getindex!(x, g::MaskedGrid, idx::Int) = getindex!(x, g.grid, g.indices[idx,:]...)

# Don't use getindex! on 1D grids
getindex{G,ID}(g::MaskedGrid{G,ID,1}, idx) = getindex(g.grid, g.indices[idx,1])



"""
An IndexedSubGrid is a subgrid corresponding to a certain range of indices of the
underlying (one-dimensional) grid.
"""
immutable IndexedSubGrid{G <: AbstractGrid1d, T} <: AbstractSubGrid{1,T}
	grid	::	G
	i1		::	Int
	i2		::	Int

	function IndexedSubGrid(grid::AbstractGrid1d{T}, i1, i2)
		@assert 1 <= i1 <= length(grid)
		@assert 1 <= i2 <= length(grid)
		@assert i1 <= i2

		new(grid, i1, i2)
	end
end

IndexedSubGrid{T}(grid::AbstractGrid1d{T}, i1, i2) = IndexedSubGrid{typeof(grid), T}(grid, i1, i2)

left(g::IndexedSubGrid) = g.grid[g.i1]

right(g::IndexedSubGrid) = g.grid[g.i2]

length(g::IndexedSubGrid) = g.i2 - g.i1 + 1

getindex(g::IndexedSubGrid, i) = g.grid[g.i1+i-1]

stepsize{G <: AbstractEquispacedGrid}(g::IndexedSubGrid{G}) = stepsize(g.grid)

range{G <: AbstractEquispacedGrid}(g::IndexedSubGrid{G}) = left(g) : stepsize(g) : right(g)





