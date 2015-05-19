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


immutable MaskedGrid{G <: AbstractGrid,ID,N,T} <: AbstractSubGrid{N,T}
	grid	::	G
	mask	::	Array{Bool,ID}
	M		::	Int

	MaskedGrid(grid::AbstractGrid{N,T}, mask) = new(grid, mask, sum(mask))
end

function MaskedGrid{N,T}(grid::AbstractGrid{N,T}, mask)
	@assert size(grid) == size(mask)

	MaskedGrid{typeof(grid),index_dim(grid),N,T}(grid, mask)
end

MaskedGrid{N,T}(grid::AbstractGrid{N,T}, domain::AbstractDomain{N,T}) = MaskedGrid(grid, in(grid, domain))

# index_dim{G,ID,N,T}(::MaskedGrid{G,ID,N,T}) = ID
# index_dim{G,ID,N,T}(::Type{MaskedGrid{G,ID,N,T}}) = ID
# index_dim{G <: MaskedGrid}(::Type{G}) = index_dim(super(G))

# TODO: add eachindex that iterates over elements that are part of the masked grid.
# eachindex_mask iterates over all elements and you have to check for element-ness manually (using `in`)
eachindex_mask(g::MaskedGrid) = eachindex(g.grid)

getindex(g::MaskedGrid, i...) = getindex(g.grid, i...)

getindex!(g::MaskedGrid, x, i...) = getindex!(g.grid, x, i...)

in(i, g::MaskedGrid) = g.mask[i]

length(g::MaskedGrid) = g.M

size(g::MaskedGrid) = (length(g),)




abstract AbstractSubIntervalGrid{G <: AbstractIntervalGrid, T} <: AbstractSubGrid{1,T}


immutable SubIntervalGrid{G <: AbstractIntervalGrid, T} <: AbstractSubIntervalGrid{G,T}
	grid	::	G
	a		::	T
	b		::	T

	function SubIntervalGrid(grid::AbstractIntervalGrid{T}, a, b)
		@assert a >= left(grid)
		@assert b <= right(grid)

		new(grid, a, b)
	end
end

SubIntervalGrid{T}(grid::AbstractIntervalGrid{T}, a, b) = SubIntervalGrid{typeof(grid), T}(grid, a, b)


immutable EquispacedSubGrid{G <: AbstractEquispacedGrid, T} <: AbstractSubIntervalGrid{G,T}
	grid	::	G
	i1		::	Int
	i2		::	Int

	function EquispacedSubGrid(grid::AbstractEquispacedGrid{T}, i1, i2)
		@assert 1 <= i1 <= length(grid)
		@assert 1 <= i2 <= length(grid)
		@assert i1 <= i2

		new(grid, i1, i2)
	end
end

EquispacedSubGrid{T}(grid::AbstractEquispacedGrid{T}, i1, i2) = EquispacedSubGrid{typeof(grid), T}(grid, i1, i2)

left(g::EquispacedSubGrid) = g.grid[g.i1]

right(g::EquispacedSubGrid) = g.grid[g.i2]

length(g::EquispacedSubGrid) = g.i2 - g.i1 + 1

stepsize(g::EquispacedSubGrid) = stepsize(g.grid)

getindex(g::EquispacedSubGrid, i) = g.grid[g.i1+i-1]

range(g::EquispacedSubGrid) = left(g) : stepsize(g) : right(g)





