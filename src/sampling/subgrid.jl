
# See also the file grid/subgrid.jl in BasisFunctions for the definition of
# AbstractSubGrid and IndexSubGrid.


"""
A MaskedGrid is a subgrid of another grid that is defined by a mask.
The mask is true or false for each point in the supergrid. The set of points
for which it is true make up the MaskedGrid.
"""
struct MaskedGrid{G,M,I,T} <: AbstractSubGrid{T,1}
    supergrid   ::	G
    mask	    ::	M
    indices     ::  Vector{I}
    M           ::	Int				# Total number of points in the mask

    MaskedGrid{G,M,I,T}(supergrid::AbstractGrid{T}, mask, indices) where {G,M,I,T} =
        new(supergrid, mask, indices, sum(mask))
end
# TODO: In MaskedGrid, perhaps we should not be storing pointers to the points of the underlying grid, but
# rather the points themselves. In that case we wouldn't need to specialize on the type of grid (parameter G can go).



function MaskedGrid(supergrid::AbstractGrid{T}, mask, indices) where {T}
	@assert size(supergrid) == size(mask)

	MaskedGrid{typeof(supergrid),typeof(mask),eltype(indices),T}(supergrid, mask, indices)
end

# These are for the assignment to indices in the function below.
convert(::Type{NTuple{N,Int}},i::CartesianIndex{N}) where {N} = ntuple(k->i[k],N)

MaskedGrid(supergrid::AbstractGrid, domain::Domain) =
	MaskedGrid(supergrid, in.(supergrid, Ref(domain)))

# MaskedGrid(maskedgrid::MaskedGrid, domain::Domain) =
#     MaskedGrid(supergrid(maskedgrid), mask(maskedgrid) .& in.(supergrid(maskedgrid), domain))

MaskedGrid(supergrid::AbstractGrid, mask) = MaskedGrid(supergrid, mask, subindices(supergrid, mask))

function subindices(supergrid, mask::BitArray)
    I = eltype(eachindex(supergrid))
    indices = Array{I}(undef, sum(mask))
    i = 1
    for m in eachindex(supergrid)
       if mask[m]
           indices[i] = m
           i += 1
       end
    end
    indices
end

size(g::MaskedGrid) = (g.M,)

mask(g::MaskedGrid) = g.mask

subindices(g::MaskedGrid) = g.indices

similar_subgrid(g::MaskedGrid, g2::AbstractGrid) = MaskedGrid(g2, g.mask, g.indices)


# Check whether element grid[i] (of the underlying grid) is in the masked grid.
issubindex(i, g::MaskedGrid) = g.mask[i]

unsafe_getindex(g::MaskedGrid, idx::Int) = unsafe_getindex(g.supergrid, g.indices[idx])

getindex(g::AbstractGrid, idx::BitArray) = MaskedGrid(g, idx)


"Create a suitable subgrid that covers a given domain."
function subgrid(grid::AbstractEquispacedGrid, domain::AbstractInterval)
    a = leftendpoint(domain)
    b = rightendpoint(domain)
    h = stepsize(grid)
    idx_a = convert(Int, ceil( (a-grid[1])/stepsize(grid))+1 )
    idx_b = convert(Int, floor( (b-grid[1])/stepsize(grid))+1 )
    idx_a = max(idx_a, 1)
    idx_b = min(idx_b, length(grid))
    IndexSubGrid(grid, idx_a:idx_b)
end

subgrid(grid::AbstractGrid, domain::Domain) = MaskedGrid(grid, domain)

function subgrid(grid::ScatteredGrid, domain::Domain)
    mask = in.(grid, Ref(domain))
    points = grid.points[mask]
    ScatteredGrid(points)
end

function subgrid(grid::MaskedGrid, domain::Domain)
    submask = in.(supergrid(grid), Ref(domain))
    MaskedGrid(supergrid(grid), submask .& mask(grid))
    # points = grid.points[mask]
    # ScatteredGrid(points)
end

subgrid(grid::ProductGrid, domain::ProductDomain) =
    ProductGrid(map(subgrid, elements(grid), elements(domain))...)

hasextension(dg::GridBasis{T,G}) where {T,G <: AbstractSubGrid} = true
hasextension(dg::GridBasis{T,G}) where {T,G <: BasisFunctions.TensorSubGrid} = true

function BasisFunctions.grid_restriction_operator(src::Dictionary, dest::Dictionary, src_grid::G, dest_grid::MaskedGrid{G,M,I,T}; options...) where {G<:AbstractGrid,M,I,T}
    @assert supergrid(dest_grid) == src_grid
    IndexRestrictionOperator(src, dest, subindices(dest_grid))
end
