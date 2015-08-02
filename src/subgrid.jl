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
    indices     ::      Array{Int,2}
    M		::	Int

    MaskedGrid(grid::AbstractGrid{N,T},  mask, indices) = new(grid, mask, indices, sum(mask))
    
end

function MaskedGrid{N,T}(grid::AbstractGrid{N,T}, mask, indices)
	@assert size(grid) == size(mask)

	MaskedGrid{typeof(grid),index_dim(grid),N,T}(grid, mask, indices)
end

function MaskedGrid{N,T}(grid::AbstractGrid{N,T}, domain::AbstractDomain{N,T})
    mask = in(grid, domain)
    indices = Array(Int,sum(mask),ndims(mask))
    i=1
    for m in eachindex(mask)
        if mask[m]
            indices[i,:]=[ind2sub(mask,m)...]
            i+=1
        end
    end
    MaskedGrid(grid, mask, indices)
end


# index_dim{G,ID,N,T}(::MaskedGrid{G,ID,N,T}) = ID
# index_dim{G,ID,N,T}(::Type{MaskedGrid{G,ID,N,T}}) = ID
# index_dim{G <: MaskedGrid}(::Type{G}) = index_dim(super(G))

length(g::MaskedGrid) = g.M

size(g::MaskedGrid) = (length(g),)


#eachindex(g::MaskedGrid) = MaskedGridRange(g, eachindex(g.grid))

# Check whether element grid[i] (of the underlying grid) is in the masked grid.
in(i, g::MaskedGrid) = g.mask[i]

immutable MaskedGridRange{G <: MaskedGrid,ITER}
	maskedgrid	::	G
	griditer	::	ITER
end

length(iter::MaskedGridRange) = length(iter.maskedgrid)

# A MaskedGridState refers to the state of the iterator of the underlying grid.
# It also holds a sequence number k, used to efficiently implement done (k == M+1).
immutable MaskedGridRangeState{T}
	gridstate	::	T		# state of the underlying grid iterator
	k			::	Int
end

start(iter::MaskedGridRange) = MaskedGridRangeState(start(iter.griditer),1)

function next(iter::MaskedGridRange, state)
	(i,gridstate) = next(iter.griditer, state.gridstate)
	while ~in(i, iter.maskedgrid)
		(i,gridstate) = next(iter.griditer, gridstate)
	end
	(i, MaskedGridRangeState(gridstate, state.k+1))
end

done(iter::MaskedGridRange, state) = (state.k == length(iter)+1)

getindex(g::MaskedGrid, idx) = getindex(g.grid, idx)

function getindex!(g::MaskedGrid, x, idx::Int)
    getindex!(g.grid, x, g.indices[idx,:]...)
end


getindex(g::MaskedGrid, idx::Int) = getindex(g.grid, g.indices[idx,:]...)


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


immutable EquispacedSubGrid{G <: AbstractGrid, T} <: AbstractGrid{1,T}
	grid	::	G
	i1		::	Int
	i2		::	Int

	function EquispacedSubGrid(grid::AbstractGrid{1,T}, i1, i2)
		@assert 1 <= i1 <= length(grid)
		@assert 1 <= i2 <= length(grid)
		@assert i1 <= i2

		new(grid, i1, i2)
	end
end

EquispacedSubGrid{T}(grid::AbstractGrid{1,T}, i1, i2) = EquispacedSubGrid{typeof(grid), T}(grid, i1, i2)

left(g::EquispacedSubGrid) = g.grid[g.i1]

right(g::EquispacedSubGrid) = g.grid[g.i2]

length(g::EquispacedSubGrid) = g.i2 - g.i1 + 1

stepsize(g::EquispacedSubGrid) = stepsize(g.grid)

getindex(g::EquispacedSubGrid, i) = g.grid[g.i1+i-1]

range(g::EquispacedSubGrid) = left(g) : stepsize(g) : right(g)





