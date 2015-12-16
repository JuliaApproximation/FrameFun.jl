# subgrid.jl


abstract AbstractSubGrid{N,T,G} <: AbstractGrid{N,T}

eltype{N,T,G}(::Type{AbstractSubGrid{N,T,G}}) = eltype(G)
eltype{G <: AbstractSubGrid}(::Type{G}) = eltype(super(G))

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



"""
An IndexedSubGrid is a subgrid corresponding to a certain range of indices of the
underlying (one-dimensional) grid.
"""
immutable IndexedSubGrid{G,T} <: AbstractSubGrid{1,T,G}
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





