# subgrid.jl

# See also the file grid/subgrid.jl in BasisFunctions for the definition of
# AbstractSubGrid and IndexSubGrid.


"""
A MaskedGrid is a subgrid of another grid that is defined by a mask.
The mask is true or false for each point in the supergrid. The set of points
for which it is true make up the MaskedGrid.
"""
struct MaskedGrid{G,M,I,T} <: AbstractSubGrid{T}
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
    MaskedGrid(supergrid, in.(supergrid, domain))

# MaskedGrid(maskedgrid::MaskedGrid, domain::Domain) =
#     MaskedGrid(supergrid(maskedgrid), mask(maskedgrid) .& in.(supergrid(maskedgrid), domain))

MaskedGrid(supergrid::AbstractGrid, mask) = MaskedGrid(supergrid, mask, subindices(supergrid, mask))

function subindices(supergrid, mask::BitArray)
    I= eltype(eachindex(supergrid))
    indices = Array{I}(sum(mask))
    i = 1
    for m in eachindex(supergrid)
       if mask[m]
           indices[i] = m
           i += 1
       end
    end
    indices
end

length(g::MaskedGrid) = g.M

size(g::MaskedGrid) = (length(g),)

mask(g::MaskedGrid) = g.mask

subindices(g::MaskedGrid) = g.indices

similar_subgrid(g::MaskedGrid, g2::AbstractGrid) = MaskedGrid(g2, g.mask, g.indices)


# Check whether element grid[i] (of the underlying grid) is in the masked grid.
is_subindex(i, g::MaskedGrid) = g.mask[i]

unsafe_getindex(g::MaskedGrid, idx) = unsafe_getindex(g.supergrid, g.indices[idx])


# Efficient extension operator
function apply!{G <: MaskedGrid}(op::Extension, dest, src::GridBasis{G}, coef_dest, coef_src)
    @assert length(coef_src) == length(src)
    @assert length(coef_dest) == length(dest)
    # @assert grid(dest) == supergrid(grid(src))

    grid1 = grid(src)
    fill!(coef_dest, 0)

    l = 0
    for i in eachindex(grid1.supergrid)
        if is_subindex(i, grid1)
            l += 1
            coef_dest[i] = coef_src[l]
        end
    end
    coef_dest
end


# Efficient restriction operator
function apply!{G <: MaskedGrid}(op::Restriction, dest::GridBasis{G}, src, coef_dest, coef_src)
    @assert length(coef_src) == length(src)
    @assert length(coef_dest) == length(dest)
    # This line below seems to allocate memory...
    # @assert grid(src) == supergrid(grid(dest))

    grid1 = grid(dest)

    l = 0
    for i in eachindex(grid1.supergrid)
        if is_subindex(i, grid1)
            l += 1
            coef_dest[l] = coef_src[i]
        end
    end
    coef_dest
end

struct NBIndexList{N}
    index::NTuple{N,Int}
    size::NTuple{N,Int}
end
NBIndexList(index::Base.CartesianIndex{N}, size) where {N}  = NBIndexList(index.I, size)
NBIndexList(index::Int, size)  = NBIndexList((index,), size)

@generated function Base.start(l::NBIndexList{N}) where {N}
	startargs = fill(-1, N)
    stopargs = fill(1, N)
	:(CartesianRange(CartesianIndex{$N}($(startargs...)), CartesianIndex{$N}($(stopargs...))), CartesianIndex{$N}($(startargs...)))
end
@generated function Base.next(l::NBIndexList{N}, state) where {N}
    t = Expr(:tuple, [:(if 1<=idx[$i]+l.index[$i]<=l.size[$i];idx[$i]+l.index[$i];elseif idx[$i]+l.index[$i]==0;l.size[$i];else; 1 ;end) for i in 1:N]...)
    return quote
        iter = state[1]
        iter_state = state[2]
        idx, iter_next_state = next(iter, iter_state)
        (sum(abs.(idx.I)) == 0) && ((idx, iter_next_state) = next(iter, iter_next_state))
        $t,(iter, iter_next_state)
    end
end
function Base.done(g::NBIndexList{N}, state) where {N}
    iter = state[1]
    iter_state = state[2]
    done(iter, iter_state)
end

"""
A Masked grid that contains the elements of grid that are on the boundary of the domain
"""
function boundary_grid(grid::AbstractGrid, domain::Domains.Domain)
    mask = boundary_mask(grid, domain);
    MaskedGrid(grid,mask);
end

boundary_grid(grid::MaskedGrid, domain::Domains.Domain) = boundary_grid(supergrid(grid), domain)


function boundary_mask(grid::AbstractGrid, domain::Domains.Domain)
    S = size(grid)
    m = BitArray(S...)#zeros(Bool, S...)
    m[:] = 0
    t = true
    for i in eachindex(grid)
        if Domains.indomain(grid[i], domain)
            t = true
            for bi in NBIndexList(i, S)
                if !Domains.indomain(grid[bi], domain)
                    t = false
                    break
                end
            end
            m[i] = !t
        end
    end
    m
end

# collect_neighbours!(mask::BitArray{N}, start_index, grid::MaskedGrid) where {N} =  collect_neighbours!(mask, start_index, copy(grid.mask))
collect_neighbours!(mask::BitArray{N}, start_index::CartesianIndex{N}, left_over::BitArray{N}) where {N} =  collect_neighbours!(mask, start_index.I, left_over)
function collect_neighbours!(mask::BitArray{N}, index::NTuple{N,Int}, left_over::BitArray{N}) where {N}
    mask[index...] = true
    left_over[index...] = false

    for i in FrameFun.NBIndexList(index, size(left_over))
        if left_over[i...]
            mask[i...] = true
            left_over[i...] = false
            collect_neighbours!(mask, i, left_over)
        end
    end
end

import Base:-, split
function -(m1::MaskedGrid, m2::MaskedGrid)
    @assert supergrid(m1)==supergrid(m2)
    MaskedGrid(supergrid(m1), m1.mask .& (.!m2.mask))
end

function split(m::MaskedGrid)

    L = length(m)
    s = supergrid(m)
    grid_mask = copy(m.mask)
    index = 0
    r = Array{AbstractGrid}(0)
    mask_total = BitArray(size(s)...)
    mask_total[:] = 0
    mask = BitArray(size(s)...)
    while index < L
        for i in subindices(m)
            if !mask_total[i]
                mask = BitArray(size(s)...)
                mask[:] = 0
                collect_neighbours!(mask, i, grid_mask)
                index += sum(mask)
                mask_total = mask_total .| mask
                push!(r, MaskedGrid(s, mask))
            end
        end
    end
    if length(r) == 1
        return [m]
    else
        return r
    end
end

function split_in_IndexGrids(m::MaskedGrid)
    L = length(m)
    s = supergrid(m)
    index = 0
    r = IndexSubGrid[]
    grid_mask = copy(m.mask)
    mask_total = BitArray(size(s)...)
    mask_total[:] = 0
    while index < L
        for i in subindices(m)
            if !mask_total[i]
                mask = BitArray(size(s)...)
                mask[:] = 0
                collect_neighbours!(mask, i, grid_mask)
                index += sum(mask)
                mask_total = mask_total .| mask

                push!(r, IndexSubGrid(s, subindices(s,mask)))
            end
        end
    end
    if length(r) == 1
        return [m]
    else
        return r
    end
end

# It is assumed that all points of `from` are in `relativeto` and that the supergrids of both grids are equal.
function relative_indices(from::MaskedGrid, relativeto::Union{IndexSubGrid,MaskedGrid})
    # @assert (supergrid(from)) == (supergrid(relativeto))
    @assert length(supergrid(from)) == length(supergrid(relativeto))
    support_index = Array{Int}(length(from))
    index = 1
    for (i_i,i) in enumerate(subindices(relativeto))
        if is_subindex(i, from)
            support_index[index] = i_i
            index += 1
        end
    end
    support_index
end

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
    mask = in.(grid, domain)
    points = grid.points[mask]
    ScatteredGrid(points)
end

function subgrid(grid::MaskedGrid, domain::Domain)
    mask = in.(supergrid(grid), domain)
    MaskedGrid(supergrid(grid), mask .& BasisFunctions.mask(grid))
    # points = grid.points[mask]
    # ScatteredGrid(points)
end

# subgrid(grid::AbstractGrid, domain::DomainBoundary) = boundary(g, domain)

# subgrid(grid::ScatteredGrid, domain::DomainBoundary) = error("Trying to take the boundary within a ScatteredGrid")


# Duck typing, v1 and v2 have to implement addition/substraction and scalar multiplication
function midpoint(v1, v2, dom::Domain, tol)
    # There has to be a midpoint
    @assert in(v2,dom) != in(v1,dom)
    if in(v2,dom)
        min=v1
        max=v2
    else
        min=v2
        max=v1
    end
    mid = NaN
    while norm(max-min) > tol
        step = (max-min)/2
        mid = min+step
        in(mid,dom) ? max=mid : min=mid
    end
    mid
end

## Avoid ambiguity (because everything >=2D is tensor but 1D is not)
function boundary{TG,T}(g::ProductGrid{TG,T},dom::Domain1d)
    println("This method being called means there is a 1D ProductGrid.")
end

function boundary{G,M}(g::MaskedGrid{G,M},dom::Domain1d)
  # TODO merge supergrid?
    boundary(grid(g),dom)
end

function boundary{TG,N,T}(g::ProductGrid{TG,T},dom::EuclideanDomain{N},tol=1e-12)
    # Initialize neighbours
    neighbours=Array{Int64}(2^N-1,N)
    # adjust columns
    for i=1:N
        # The very first is not a neighbour but the point itself.
        for j=2:2^N
            neighbours[j-1,i]=(floor(Int,(j-1)/(2^(i-1))) % 2)
        end
    end
    CartesianNeighbours = Array{CartesianIndex{N}}(2^N-1)
    for j=1:2^N-1
        CartesianNeighbours[j]=CartesianIndex{N}(neighbours[j,:]...)
    end
    midpoints = eltype(g)[]
    # for each element
    for i in eachindex(g)
        # for all neighbours
        for j=1:2^N-1
            neighbour = i + CartesianNeighbours[j]
            # check if any are on the other side of the boundary
            try
                if in(g[i],dom) != in(g[neighbour],dom)
                    # add the midpoint to the grid
                    push!(midpoints, midpoint(g[i],g[neighbour],dom,tol))
                end
            catch y
                isa(y,BoundsError) || rethrow(y)
            end
        end
    end
    ScatteredGrid(midpoints)
end


function boundary(g::AbstractGrid{T},dom::Domain1d,tol=1e-12) where {T <: Number}
    midpoints = T[]
    # for each element
    for i in eachindex(g)
        # check if any are on the other side of the boundary
        try
            if in(g[i],dom) != in(g[i+1],dom)
                # add the midpoint to the grid
                push!(midpoints, midpoint(g[i],g[i+1],dom,tol))
            end
        catch y
            isa(y,BoundsError) || rethrow(y)
        end
    end
    ScatteredGrid(midpoints)
end


function boundary{G,M,N}(g::MaskedGrid{G,M},dom::EuclideanDomain{N})
    boundary(supergrid(g),dom)
end

has_extension{G <: AbstractSubGrid}(dg::DiscreteGridSpace{G}) = true

function BasisFunctions.grid_restriction_operator(src::Span, dest::Span, src_grid::G, dest_grid::MaskedGrid{G,M,I,T}; options...) where {G<:AbstractGrid,M,I,T}
    @assert supergrid(dest_grid) == src_grid
    IndexRestrictionOperator(src, dest, subindices(dest_grid))
end
