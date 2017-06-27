# extensions.jl
# A collection of extensions to the Domains package.

###########################
# Applying broadcast to in
###########################

# Intercept a broadcasted call to indomain. We assume that the user wants evaluation
# in a set of points (which we call a grid), rather than in a single point.
broadcast(::typeof(in), grid, d::Domain) = indomain_broadcast(grid, d)

# # Default methods for evaluation on a grid: the default is to call eval on the domain with
# # points as arguments. Domains that have faster grid evaluation routines may define their own version.
indomain_broadcast(grid, d::Domain) = indomain_broadcast!(zeros(Bool, size(grid)), grid, d)
# TODO: use BitArray here

function indomain_broadcast!(result, grid, domain::Domain)
    for (i,x) in enumerate(grid)
        result[i] = indomain(x, domain)
    end
    result
end

function indomain_broadcast(grid, d::UnionDomain)
    z = indomain_broadcast(grid, element(d,1))
    for i in 2:nb_elements(d)
        z = z .| indomain_broadcast(grid, element(d,2))
    end
    z
end

function indomain_broadcast(grid, d::IntersectionDomain)
    z = indomain_broadcast(grid, element(d,1))
    for i in 2:nb_elements(d)
        z = z .& indomain_broadcast(grid, element(d,2))
    end
    z
end

function indomain_broadcast(grid, d::DifferenceDomain)
    z1 = indomain_broadcast(grid, d.d1)
    z2 = indomain_broadcast(grid, d.d2)
    z1 .& (.~z2)
end

# This breaks mappeddomain in for grids

## indomain_broadcast(grid, d::DerivedDomain) = indomain_broadcast(grid, superdomain(d))
## indomain_broadcast!(result, grid, d::DerivedDomain) = indomain_broadcast!(result, grid, superdomain(d))

# Check whether a value is in an interval, up to 10 times machine precision
in(x::Number, a::T, b::T) where {T <: AbstractFloat} = (a-10eps(T) <= x <= b+10eps(T))
in(x::Number, a::T, b::T) where {T <: Number} = a <= x <= b


#################
# Bounding boxes
#################

boundingbox(d::AbstractInterval) = BBox(leftendpoint(d), rightendpoint(d))

boundingbox(::UnitBall{N,T}) where {N,T} = BBox{N,T}(-ones(SVector{N,T}), ones(SVector{N,T}))

boundingbox(c::Ball) = BBox((c.center[1]-c.radius,c.center[2]-c.radius,c.center[3]-c.radius),(c.center[1]+c.radius,c.center[2]+c.radius,c.center[3]+c.radius))

boundingbox(d::ProductDomain) = tensorproduct(map(boundingbox, elements(d))...)

boundingbox(d::DerivedDomain) = boundingbox(superdomain(d))

boundingbox(d::UnionDomain) = ∪(map(boundingbox, elements(d))...)

boundingbox(d::IntersectionDomain) = ∩(map(boundingbox, elements(d))...)

boundingbox(d::DifferenceDomain) = boundingbox(d.d1)




# Now here is a problem: how do we compute a bounding box, without extra knowledge
# of the map? We can only do this for some maps.
boundingbox(d::MappedDomain) = mapped_boundingbox(boundingbox(superdomain(d)), mapping(d))

function mapped_boundingbox(box::BBox1, fmap)
    l,r = box[1]
    ml = fmap*l
    mr = fmap*r
    BBox(min(ml,mr), max(ml,mr))
end

# In general, we can at least map all the corners of the bounding box of the
# underlying domain, and compute a bounding box for those points. This will be
# correct for affine maps.
function mapped_boundingbox(box::BBox{N}, fmap) where {N}
    crn = corners(box)
    mapped_corners = [fmap*c for c in crn]
    left = [minimum([mapped_corners[i][j] for i in 1:length(mapped_corners)]) for j in 1:N]
    right = [maximum([mapped_corners[i][j] for i in 1:length(mapped_corners)]) for j in 1:N]
    BBox(left, right)
end

# We can do better for diagonal maps, since the problem simplifies: each dimension
# is mapped independently.
mapped_boundingbox(box::BBox{N}, fmap::ProductMap) where {N} =
    tensorproduct([mapped_boundingbox(element(box,i), element(fmap,i)) for i in 1:N]...)

boundingbox(d::TranslatedDomain) = boundingbox(domain(d)) + translationvector(d)

##########################################################################
### Distances and Normals
##########################################################################

# Function that returns the normal vector when evaluated at the domain boundary
normal(x, d::Domain) = error("Normal not available for this domain type")

# Auxiliary function that returns distance from the boundary in some metric
dist(x, d::Domain) = error("Domain distance not available for this domain type")

dist(x, ::UnitSimplex) = min(minimum(x),1-sum(x))

function normal(x, ::UnitSimplex)
    if minimum(x)<abs(sum(x)-1)/sqrt(length(x))
        index = findmin(x)[2]
        z = zeros(x)
        setindex(z,-1,index)
    else
        ones(x)
    end
end
normal(x, d::UnitBall) = x/norm(x)

dist(x, d::AbstractInterval) = min(rightendpoint(d)-x,x-leftendpoint(d))

normal(x, d::AbstractInterval) = abs(leftendpoint(d)-x) < abs(rightendpoint(d)-x) ? -1:1

dist(x, d::UnitBall) = 1-norm(x)

normal(x, d::UnitSphere) = x/norm(x)

dist(x, d::UnitSphere) = 1-norm(x)

dist(x,d::UnionDomain) = indomain(x,d) ? sum(map(di->max(0,dist(x,di)),elements(d))) : maximum(map(di->dist(x,di),elements(d)))

normal(x,d::UnionDomain) = normal(x,elements(d)[findmin(map(di->abs(dist(x,di)),elements(d)))[2]])

dist(x,d::IntersectionDomain) = minimum(map(di->dist(x,di),elements(d)))

normal(x,d::IntersectionDomain) = normal(x,elements(d)[findmin(map(di->dist(x,di),elements(d)))[2]])

dist(x,d::DifferenceDomain) = indomain(x,d) ? min(abs(dist(x,d.d1)),abs(dist(x,d.d2))) : -1*min(abs(dist(x,d.d1)),abs(dist(x,d.d2)))

normal(x,d::DifferenceDomain) = dist(x,d.d1)<dist(x,d.d2) ? normal(x,d.d1) : -1*normal(x,d.d2)

dist(x, t::ProductDomain) = minimum(map(dist,x,elements(t)))

dist(x, d::MappedDomain) = dist(mapping(d)*x,superdomain(d))

function normal(x, d::MappedDomain)
    x = applymap(mapping(d),normal(mapping(d)*x,superdomain(d)))
    x0 = apply_inverse(mapping(d),zeros(size(x)))
   (x-x0)/norm(x-x0)
end
function normal(x, t::ProductDomain)
    index = findmin(map(dist,x,elements(t)))[2]
    [(i==index)*normal(x[i],element(t,i)) for i =1:length(elements(t))]
end
##########################################################################
### Assorted Domains
##########################################################################


function randomcircles(n, radius = 0.3)
    list = [disk(radius, SVector(((2*rand(2)-1)*0.8)...)) for i=1:n]
    UnionDomain(list...)
end

ndims(::Type{Domain{T}}) where {T} = ndims_type(T)
ndims(::Type{D}) where {D <: Domain} = ndims(supertype(D))
ndims(d::Domain) = ndims(typeof(d))
ndims_type(::Type{SVector{N,T}}) where {N,T} = N
ndims_type(::Type{T}) where {T <: Number} = 1
