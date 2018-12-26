# A collection of extensions to the DomainSets package.

using DomainSets: inverse_map, forward_map

###########################
# Applying broadcast to in
###########################

# Intercept a broadcasted call to indomain. We assume that the user wants evaluation
# in a set of points (which we call a grid), rather than in a single point.
# TODO: the user may want to evaluate a single point in a sequence of domains...
broadcast(::typeof(in), grid, d::Domain) = indomain_broadcast(grid, d)

# # Default methods for evaluation on a grid: the default is to call eval on the domain with
# # points as arguments. Domains that have faster grid evaluation routines may define their own version.
indomain_broadcast(grid, d::Domain) = indomain_broadcast!(BitArray(undef, size(grid)), grid, d)
# TODO: use BitArray here

function indomain_broadcast!(result, grid, domain::Domain)
    for (i,x) in enumerate(grid)
        result[i] = DomainSets.indomain(x, domain)
    end
    result
end

function indomain_broadcast(grid, d::UnionDomain)
    z = indomain_broadcast(grid, element(d,1))
    for i in 2:numelements(d)
        z = z .| indomain_broadcast(grid, element(d,i))
    end
    z
end

function indomain_broadcast(grid, d::IntersectionDomain)
    z = indomain_broadcast(grid, element(d,1))
    for i in 2:numelements(d)
        z = z .& indomain_broadcast(grid, element(d,i))
    end
    z
end

function indomain_broadcast(grid, d::DifferenceDomain)
    z1 = indomain_broadcast(grid, d.d1)
    z2 = indomain_broadcast(grid, d.d2)
    z1 .& (.~z2)
end

# This breaks mappeddomain in for grids

## indomain_broadcast(grid, d::DerivedDomain) = indomain_broadcast(grid, source(d))
## indomain_broadcast!(result, grid, d::DerivedDomain) = indomain_broadcast!(result, grid, souce(d))

# Check whether a value is in an interval, up to 10 times machine precision
in(x::Number, a::T, b::T) where {T <: AbstractFloat} = (a-10eps(T) <= x <= b+10eps(T))
in(x::Number, a::T, b::T) where {T <: Number} = a <= x <= b


#################
# Bounding boxes
#################

# A bounding box is an Interval or ProductDomain of intervals that encompasses the domain.

# If the boundingbox is not a product of intervals, something has gone wrong.

# Some of these constructors can hopefully disappear when spaces are introduced.

boundingbox(a::SVector{1}, b::SVector{1}) = a[1]..b[1]

boundingbox(a::Number, b::Number) = a..b

boundingbox(a, b) = cube(a,b)

boundingbox(d::AbstractInterval) = d

boundingbox(::UnitHyperBall{N,T}) where {N,T} = cube(-ones(SVector{N,T}), ones(SVector{N,T}))

boundingbox(d::ProductDomain) = cartesianproduct(map(boundingbox, elements(d))...)

boundingbox(d::DerivedDomain) = boundingbox(source(d))

boundingbox(d::DifferenceDomain) = boundingbox(d.d1)

# Extra functions

minimum(d::Domain) = minimum(boundingbox(d))

maximum(d::Domain) = maximum(boundingbox(d))

minimum(box::ProductDomain) = SVector(map(minimum,elements(box)))

maximum(box::ProductDomain) = SVector(map(maximum,elements(box)))

# TODO: improve when grids move into domains?
equispaced_grid(d::Domain, ns) = cartesianproduct([PeriodicEquispacedGrid(ns[idx], infimum(d)[idx], supremum(d)[idx]) for idx = 1:dimension(boundingbox(d))]...)

function boundingbox(d::UnionDomain)
    left = SVector(minimum(hcat(map(infimum,elements(d))...);dims=2)...)
    right = SVector(maximum(hcat(map(supremum,elements(d))...);dims=2)...)
    boundingbox(left,right)
end

function boundingbox(d::IntersectionDomain)
    left = SVector(maximum(hcat(map(infimum,elements(d))...);dims=2)...)
    right = SVector(minimum(hcat(map(supremum,elements(d))...);dims=2)...)
    boundingbox(left,right)
end

DomainSets.superdomain(d::DomainSets.MappedDomain) = DomainSets.source(d)

# Now here is a problem: how do we compute a bounding box, without extra knowledge
# of the map? We can only do this for some maps.
boundingbox(d::DomainSets.MappedDomain) = mapped_boundingbox(boundingbox(source(d)), forward_map(d))

function mapped_boundingbox(box::Interval, fmap)
    l,r = (minimum(box),maximum(box))
    ml = fmap*l
    mr = fmap*r
    boundingbox(min(ml,mr), max(ml,mr))
end

# In general, we can at least map all the corners of the bounding box of the
# underlying domain, and compute a bounding box for those points. This will be
# correct for affine maps.
function mapped_boundingbox(box::ProductDomain, fmap)
    crn = corners(minimum(box),maximum(box))
    mapped_corners = [fmap*crn[:,i] for i in 1:size(crn,2)]
    left = [minimum([mapped_corners[i][j] for i in 1:length(mapped_corners)]) for j in 1:size(crn,1)]
    right = [maximum([mapped_corners[i][j] for i in 1:length(mapped_corners)]) for j in 1:size(crn,1)]
    cube(left, right)
end

# Auxiliary functions to rotate a bounding box when mapping it.
function corners(left::AbstractVector, right::AbstractVector)
    @assert length(left)==length(right)
    N=length(left)
    corners = zeros(N,2^N)
    # All possible permutations of the corners
    for i=1:2^length(left)
        for j=1:N
            corners[j,i] = ((i>>(j-1))%2==0) ? left[j] : right[j]
        end
    end
    corners
end


##########################################################################
### Distances and Normals
##########################################################################

# Function that returns the normal vector when evaluated at the domain boundary
# This vector is of unit length and points to the OUTSIDE.
normal(x, d::Domain) = error("Normal not available for this domain type")

# Auxiliary function that returns distance from the boundary in some metric
distance(x, d::Domain) = error("Domain distance not available for this domain type")

distance(x, ::UnitSimplex) = min(minimum(x),1-sum(x))

function normal(x, ::UnitSimplex)
    z = fill(eltype(x)(0), size(x))
    if minimum(x)<abs(sum(x)-1)/sqrt(length(x))
        index = findmin(x)[2]
        setindex!(z,-1,index)
    else
        z = fill(eltype(x)(1), size(x))
    end
    return z/norm(z)
end

normal(x, d::UnitHyperBall) = x/norm(x)

distance(x, d::AbstractInterval) = min(maximum(d)-x,x-minimum(d))

normal(x, d::AbstractInterval) = (abs(minimum(d)-x) < abs(maximum(d)-x)) ? -1 : 1

distance(x, d::UnitHyperBall) = 1-norm(x)

normal(x, d::UnitHyperSphere) = x/norm(x)

distance(x, d::UnitHyperSphere) = 1-norm(x)

distance(x,d::UnionDomain) = indomain(x,d) ? sum(map(di->max(0,distance(x,di)),elements(d))) : maximum(map(di->distance(x,di),elements(d)))

normal(x,d::UnionDomain) = normal(x,elements(d)[findmin(map(di->abs(distance(x,di)),elements(d)))[2]])

distance(x,d::IntersectionDomain) = minimum(map(di->distance(x,di),elements(d)))

normal(x,d::IntersectionDomain) = normal(x,elements(d)[findmin(map(di->distance(x,di),elements(d)))[2]])

distance(x,d::DifferenceDomain) = indomain(x,d) ? min(abs(distance(x,d.d1)),abs(distance(x,d.d2))) : -1*min(abs(distance(x,d.d1)),abs(distance(x,d.d2)))

normal(x,d::DifferenceDomain) = abs(distance(x,d.d1))<abs(distance(x,d.d2)) ? normal(x,d.d1) : -1*normal(x,d.d2)

distance(x, t::ProductDomain) = minimum(map(distance, x, elements(t)))

distance(x, d::DomainSets.MappedDomain) = distance(inverse_map(d)*x,source(d))

function normal(x, d::DomainSets.MappedDomain)
    x = applymap(inverse_map(d),normal(inverse_map(d)*x,source(d)))
    x0 = apply_inverse(inverse_map(d),zeros(size(x)))
   (x-x0)/norm(x-x0)
end
function normal(x, t::ProductDomain)
    index = findmin(map(distance, x, elements(t)))[2]
    [(i==index)*normal(x[i],element(t,i)) for i =1:length(elements(t))]
end


##########################################################################
### Volume
##########################################################################

volume(d::Domain) = missing

volume(d::AbstractInterval) = abs(rightendpoint(d)-leftendpoint(d))
