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
function boundary(g::ProductGrid{TG,T,N},dom::Domain1d) where {TG,T,N}
    println("This method being called means there is a 1D ProductGrid.")
end

function boundary(g::MaskedGrid{G,M},dom::Domain1d) where {G,M}
  # TODO merge supergrid?
    boundary(grid(g),dom)
end

function boundary(g::ProductGrid{TG,T},dom::EuclideanDomain{N},tol=1e-12) where {TG,N,T}
    # Initialize neighbours
    neighbours = Array{Int64}(undef, 2^N-1,N)
    # adjust columns
    for i=1:N
        # The very first is not a neighbour but the point itself.
        for j=2:2^N
            neighbours[j-1,i]=(floor(Int,(j-1)/(2^(i-1))) % 2)
        end
    end
    CartesianNeighbours = Array{CartesianIndex{N}}(undef,2^N-1)
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


function boundary(g::MaskedGrid{G,M},dom::EuclideanDomain{N}) where {G,M,N}
    boundary(supergrid(g),dom)
end
