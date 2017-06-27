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

indomain_broadcast(grid, d::DerivedDomain) = indomain_broadcast(grid, superdomain(d))
indomain_broadcast!(result, grid, d::DerivedDomain) = indomain_broadcast!(result, grid, superdomain(d))

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
