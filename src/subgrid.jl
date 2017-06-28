# subgrid.jl

# See also the file grid/subgrid.jl in BasisFunctions for the definition of
# AbstractSubGrid and IndexSubGrid.

"""
A MaskedGrid is a subgrid of another grid that is defined by a mask.
The mask is true or false for each point in the supergrid. The set of points
for which it is true make up the MaskedGrid.
"""
struct MaskedGrid{G,M,I,N,T} <: AbstractSubGrid{N,T}
    supergrid   ::	G
    mask	    ::	M
    indices     ::  Vector{I}
    M           ::	Int				# Total number of points in the mask

    MaskedGrid{G,M,I,N,T}(supergrid::AbstractGrid{N,T}, mask, indices) where {G,M,I,N,T} =
        new(supergrid, mask, indices, sum(mask))
end
# TODO: In MaskedGrid, perhaps we should not be storing pointers to the points of the underlying grid, but
# rather the points themselves. In that case we wouldn't need to specialize on the type of grid (parameter G can go).

function MaskedGrid{N,T}(supergrid::AbstractGrid{N,T}, mask, indices)
	@assert size(supergrid) == size(mask)

	MaskedGrid{typeof(supergrid),typeof(mask),eltype(indices),N,T}(supergrid, mask, indices)
end

# These are for the assignment to indices in the function below.
convert{N}(::Type{NTuple{N,Int}},i::CartesianIndex{N}) = ntuple(k->i[k],N)

function MaskedGrid(supergrid::AbstractGrid, domain::Domain)
    mask = in.(supergrid, domain)
    I = eltype(eachindex(supergrid))
    indices = Array{I}(sum(mask))
    i = 1
    for m in eachindex(supergrid)
        if mask[m]
            indices[i] = m
            i += 1
        end
    end
    MaskedGrid(supergrid, mask, indices)
end


length(g::MaskedGrid) = g.M

size(g::MaskedGrid) = (length(g),)

similar_subgrid(g::MaskedGrid, g2::AbstractGrid) = MaskedGrid(g2, g.mask, g.indices)


# Check whether element grid[i] (of the underlying grid) is in the masked grid.
is_subindex(i, g::MaskedGrid) = g.mask[i]

unsafe_getindex(g::MaskedGrid, idx) = unsafe_getindex(g.supergrid, g.indices[idx])


# Efficient extension operator
function apply!{G <: MaskedGrid}(op::Extension, dest, src::DiscreteGridSpace{G}, coef_dest, coef_src)
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
function apply!{G <: MaskedGrid}(op::Restriction, dest::DiscreteGridSpace{G}, src, coef_dest, coef_src)
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
function boundary{TG,T}(g::TensorProductGrid{TG,1,T},dom::Domain1d)
    println("This method being called means there is a 1D tensorproductgrid.")
end

function boundary{G,M}(g::MaskedGrid{G,M,1},dom::Domain1d)
  # TODO merge supergrid?
    boundary(grid(g),dom)
end

function boundary{TG,N,T}(g::TensorProductGrid{TG,N,T},dom::EuclideanDomain{N},tol=1e-12)
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
    midpoints = SVector{N,T}[]
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


function boundary{T}(g::AbstractGrid{1,T},dom::Domain1d,tol=1e-12)
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


function boundary{G,M,N}(g::MaskedGrid{G,M,N},dom::EuclideanDomain{N})
    boundary(supergrid(g),dom)
end

# function evaluation_operator{G <: AbstractSubGrid}(s::FunctionSet, d::DiscreteGridSpace{G})
#     d2 = DiscreteGridSpace(grid(grid(d)), eltype(s))
#     restriction_operator(d2, d) * evaluation_operator(s, d2)
# end

has_extension{G <: AbstractSubGrid}(dg::DiscreteGridSpace{G}) = true
