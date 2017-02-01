# derived_domains.jl


################################################################################
### A domain described by a characteristic function
################################################################################

immutable Characteristic{N,T} <: AbstractDomain{N}
    char    ::  Function
    box    ::  BBox{N,T}
end

Characteristic{N}(char::Function, dom::AbstractDomain{N}) = Characteristic(char,boundingbox(dom))

indomain(x, c::Characteristic) = c.char(x)

boundingbox(c::Characteristic) = c.box

show(io::IO, c::Characteristic) = print(io, "a domain described by a characteristic function")


################################################################################
### A domain described by a grid
################################################################################

# Implicit assumption: the supergrid of G is an equispaced 1d grid, or a tensorproduct of those grids.
immutable EquispacedGridDomain{N,T} <: AbstractDomain{N} 
    G    ::  AbstractGrid
end

EquispacedGridDomain{N,T}(g::AbstractGrid{N,T}) = EquispacedGridDomain{N,T}(g)

# For indomain we assume the grid is tensorproduct if multidimensional (only to be used with function values from arrays)
function indomain(x, gd::EquispacedGridDomain)
    index = (x-left(supergrid(gd.G)))./(right(supergrid(gd.G))-left(supergrid(gd.G)))
    N = ndims(gd)
    neighbours=Array(Int64,N)
    # Check the bounding rectangle
    for j=1:2^N
        for i=1:1
            neighbours[i]=(floor(Int,(j-1)/(2^(i-1))) % 2)
        end
        mdindex = ceil(Int,index.*(SVector(size(supergrid(gd.G)))))+neighbours
        # Check for bounding rectangle 
        if (prod(mdindex.>0) && prod(mdindex.<=[size(supergrid(gd.G))...])) 
            rindex = linear_index(supergrid(gd.G),Tuple(mdindex))
            is_subindex(rindex,gd.G) && return true
        end
    end
    # All bounding rectangle points contain data
    return false
end

boundingbox(gd::EquispacedGridDomain) = BBox(left(supergrid(gd.G)),right(supergrid(gd.G)))

show(io::IO, gd::EquispacedGridDomain) = print(io, "a domain described by a subgrid of some equispaced grid")



################################################################################
### The union of two domains
################################################################################

immutable DomainUnion{D1,D2,N} <: AbstractDomain{N}
    d1    ::  D1
    d2    ::  D2

    DomainUnion(d1::AbstractDomain{N}, d2::AbstractDomain{N}) = new(d1, d2)
end

DomainUnion{N}(d1::AbstractDomain{N}, d2::AbstractDomain{N}) = DomainUnion{typeof(d1),typeof(d2),N}(d1, d2)

union(d1::AbstractDomain, d2::AbstractDomain) = (d1 == d2 ? d1 : DomainUnion(d1,d2))


# The union of two domains corresponds to a logical OR of their characteristic functions
indomain(x, d::DomainUnion) = in(x, d.d1) || in(x, d.d2)

function indomain_grid(g::AbstractGrid, d::DomainUnion)
    z1 = indomain_grid(g, d.d1)
    z2 = indomain_grid(g, d.d2)
    z1 | z2
end

(+)(d1::AbstractDomain, d2::AbstractDomain) = union(d1,d2)
(|)(d1::AbstractDomain, d2::AbstractDomain) = union(d1,d2)


boundingbox(d::DomainUnion) = boundingbox(d.d1) ∪ boundingbox(d.d2)

function show(io::IO, d::DomainUnion)
    print(io, "a union of two domains: \n")
    print(io, "    First domain: ", d.d1, "\n")
    print(io, "    Second domain: ", d.d2, "\n")
end


################################################################################
### The intersection of two domains
################################################################################

immutable DomainIntersection{D1,D2,N} <: AbstractDomain{N}
    d1    ::  D1
    d2    ::  D2

    DomainIntersection(d1::AbstractDomain{N}, d2::AbstractDomain{N}) = new(d1, d2)
end

DomainIntersection{N}(d1::AbstractDomain{N},d2::AbstractDomain{N}) = DomainIntersection{typeof(d1),typeof(d2),N}(d1, d2)

# The intersection of two domains corresponds to a logical AND of their characteristic functions
indomain(x, d::DomainIntersection) = in(x, d.d1) && in(x, d.d2)

function indomain_grid(g::AbstractGrid, d::DomainIntersection)
    z1 = indomain_grid(g, d.d1)
    z2 = indomain_grid(g, d.d2)
    z1 & z2
end

(&)(d1::AbstractDomain, d2::AbstractDomain) = intersect(d1,d2)

intersect(d1::AbstractDomain, d2::AbstractDomain) = (d1 == d2 ? d1 : DomainIntersection(d1,d2))

function intersect(d1::TensorProductDomain, d2::TensorProductDomain)
    @assert ndims(d1) == ndims(d2)
    if composite_length(d1) == composite_length(d2)
        tensorproduct([intersect(element(d1,i), element(d2,i)) for i in 1:composite_length(d1)]...)
    else
        DomainIntersection(d1, d2)
    end
end



boundingbox(d::DomainIntersection) = boundingbox(d.d1) ∩ boundingbox(d.d2)

function show(io::IO, d::DomainIntersection)
    print(io, "the intersection of two domains: \n")
    print(io, "    First domain: ", d.d1, "\n")
    print(io, "    Second domain: ", d.d1, "\n")
end


################################################################################
### The difference between two domains
################################################################################

immutable DomainDifference{D1,D2,N} <: AbstractDomain{N}
    d1    ::  D1
    d2    ::  D2

    DomainDifference(d1::AbstractDomain{N}, d2::AbstractDomain{N}) = new(d1, d2)
end

DomainDifference{N}(d1::AbstractDomain{N}, d2::AbstractDomain{N}) = DomainDifference{typeof(d1),typeof(d2),N}(d1,d2)

setdiff(d1::AbstractDomain, d2::AbstractDomain) = DomainDifference(d1, d2)

# The difference between two domains corresponds to a logical AND NOT of their characteristic functions
indomain(x, d::DomainDifference) = indomain(x, d.d1) && (~indomain(x, d.d2))

function indomain_grid(g::AbstractGrid, d::DomainDifference)
    z1 = indomain_grid(g, d.d1)
    z2 = indomain_grid(g, d.d2)
    z1 & (~z2)
end

(-)(d1::AbstractDomain, d2::AbstractDomain) = setdiff(d1, d2)
(\ )(d1::AbstractDomain, d2::AbstractDomain) = setdiff(d1, d2)


boundingbox(d::DomainDifference) = boundingbox(d.d1)

function show(io::IO, d::DomainDifference)
    print(io, "the difference of two domains: \n")
    print(io, "    First domain: ", d.d1, "\n")
    print(io, "    Second domain: ", d.d2, "\n")
end


################################################################################
### A revolved domain is a 2D-domain rotated about the X-axis
################################################################################

immutable RevolvedDomain{D} <: AbstractDomain{3}
    d     ::  D
end

revolve(d::AbstractDomain{2}) = RevolvedDomain(d)

function indomain(x, d::RevolvedDomain)
    r = sqrt(x[2]^2+x[3])
    phi = atan2(x[2]/x[1])
    theta = acos(x[3]/r)
    indomain((x[1],r), d.d)
end


boundingbox(d::RevolvedDomain) = BBox((left(d.d)[1],left(d.d)...),(right(d.d)[1],right(d.d)...))

function show(io::IO, r::RevolvedDomain)
    print(io, "the revolution of: ", r.d1)
end


################################################################################
### A rotated domain
################################################################################

immutable RotatedDomain{D,T,N,L} <: AbstractDomain{N}
   d                 ::  D
   angle             ::  Vector{T}
   rotationmatrix    ::  SMatrix{N,N,T,L}

    RotatedDomain(d,angle,rotationmatrix) = new(d, angle, rotationmatrix)
end

# Rotation in positive direction
rotationmatrix(theta) = SMatrix{2,2}([cos(theta) -sin(theta); sin(theta) cos(theta)])
# Rotation about X-axis (phi), Y-axis (theta) and Z-axis (psi)
rotationmatrix(phi,theta,psi) =
   SMatrix{3,3}([cos(theta)*cos(psi) cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi) sin(phi)*sin(psi)-cos(phi)*sin(theta)*cos(psi); -cos(theta)*sin(psi) cos(phi)*cos(psi)-sin(phi)*sin(theta)*sin(psi) sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi); sin(theta) -sin(phi)*cos(theta) cos(phi)*cos(theta)])

RotatedDomain{N,T}(d::AbstractDomain{N}, angle::Vector{T}, m::SMatrix{N,N,T} = rotationmatrix(theta)) =
   RotatedDomain{typeof(d),T,N,N*N}(d, angle, m)

RotatedDomain(d::AbstractDomain{2}, theta::Number) = RotatedDomain{2,typeof(theta),typeof(d)}(d, [theta], rotationmatrix(theta))
# types annotated to remove ambiguity
RotatedDomain{T,D}(d::D, phi::T, theta::T, psi::T) = RotatedDomain{3,T,D}(d, [phi,theta,psi], rotationmatrix(phi,theta,psi))

rotate{T}(d::AbstractDomain{2}, theta::T) = RotatedDomain(d, theta)
rotate{T}(d::AbstractDomain{3}, phi::T, theta::T, psi::T) = RotatedDomain(d, phi, theta, psi)

indomain(x, d::RotatedDomain) = indomain(d.rotationmatrix*x, d.d)

(==)(d1::RotatedDomain, d2::RotatedDomain) = (d1.d == d2.d) && (d1.angle == d2.angle) #&& (d1.rotationmatrix == d2.rotationmatrix)

# very crude bounding box (doesn't work!!!)
boundingbox(r::RotatedDomain)= sqrt(2)*boundingbox(r.d)




################################################################################
### A translated domain
################################################################################

immutable TranslatedDomain{D,T,N} <: AbstractDomain{N}
    domain  ::  D
    trans   ::  SVector{N,T}

    TranslatedDomain(domain::AbstractDomain{N}, trans) = new(domain, trans)
end

TranslatedDomain{N}(domain::AbstractDomain{N}, trans::SVector{N}) = TranslatedDomain{typeof(domain),eltype(trans),N}(domain, trans)

domain(d::TranslatedDomain) = d.domain

translationvector(d::TranslatedDomain) = d.trans

function indomain(x, d::TranslatedDomain)
    indomain(x-d.trans, d.domain)
end

(+)(d::AbstractDomain, trans::SVector) = TranslatedDomain(d, trans)
(+)(d::TranslatedDomain, trans::SVector) = TranslatedDomain(domain(d), trans+translationvector(d))

boundingbox(d::TranslatedDomain) = boundingbox(domain(d)) + translationvector(d)



################################################################################
### A collection of domains
################################################################################

type DomainCollection{N} <: AbstractDomain{N}
    list    ::  Array{AbstractDomain{N},1}
end

DomainCollection(d::AbstractDomain) = DomainCollection([d])


length(d::DomainCollection) = length(d.list)

domain(d::DomainCollection, i) = d.list[i]

# Iteration over the domain list
start(d::DomainCollection) = start(d.list)

next(d::DomainCollection, state) = next(d.list, state)

done(d::DomainCollection, state) = done(d.list, state)

function indomain(x, dc::DomainCollection)
    z = false
    for d in dc
        z = z || indomain(x, d)
    end
    z
end

function indomain_grid!(z, g::AbstractGrid, dc::DomainCollection)
    for d in dc
        indomain_grid!(z, g, d)
    end
    z
end

push!(dc::DomainCollection, d::AbstractDomain) = push!(dc.list, d)



function boundingbox(d::DomainCollection)
    ubox = boundingbox(d.list[1])
    for i = 2:length(d.list)
        ubox = union(ubox, boundingbox(d.list[i]))
    end
    ubox
end


show(io::IO, d::DomainCollection) = print(io, "a collection of ", length(d.list), " domains")


#######################
# Mapped domains
#######################

"""
A MappedDomain consists of a domain and a bidirectional map. The forward map
maps the domain onto the mapped domain, and the inverse map maps it back.
A point lies in the mapped domain, if the inverse map of that point lies in the
original domain.
"""
# TODO: experiment with leaving out the type parameters and implement fast indomain_grid
immutable MappedDomain{DOMAIN <: AbstractDomain,MAP,N} <: AbstractDomain{N}
    domain  ::  DOMAIN
    # The forward map, from the underlying domain to the mapped domain
    fmap    ::  MAP

    # With this inner constructor we enforce that N is the dimension of the domain
    MappedDomain(domain::AbstractDomain{N}, fmap) = new(domain, fmap)
end

MappedDomain{N}(domain::AbstractDomain{N}, fmap) =
    MappedDomain{typeof(domain),typeof(fmap),N}(domain, fmap)

domain(d::MappedDomain) = d.domain

mapping(d::MappedDomain) = d.fmap

indomain(x, d::MappedDomain) = indomain(inverse_map(mapping(d), x), domain(d))

# Now here is a problem: how do we compute a bounding box, without extra knowledge
# of the map? We can only do this for some maps.
boundingbox(d::MappedDomain) = mapped_boundingbox(boundingbox(domain(d)), mapping(d))

function mapped_boundingbox(box::BBox1, fmap)
    l,r = box[1]
    ml = fmap*l
    mr = fmap*r
    BBox(min(ml,mr), max(ml,mr))
end

# In general, we can at least map all the corners of the bounding box of the
# underlying domain, and compute a bounding box for those points. This will be
# correct for affine maps.
function mapped_boundingbox{N}(box::BBox{N}, fmap)
    crn = corners(box)
    mapped_corners = [fmap*c for c in crn]
    left = [minimum([mapped_corners[i][j] for i in 1:length(mapped_corners)]) for j in 1:N]
    right = [maximum([mapped_corners[i][j] for i in 1:length(mapped_corners)]) for j in 1:N]
    BBox(left, right)
end

# We can do better for diagonal maps, since the problem simplifies: each dimension
# is mapped independently.
mapped_boundingbox{N}(box::BBox{N}, fmap::DiagonalMap) =
    tensorproduct([mapped_boundingbox(element(box,i), element(fmap,i)) for i in 1:N]...)

apply_map(domain::AbstractDomain, map::AbstractMap) = MappedDomain(domain, map)

apply_map(d::MappedDomain, map::AbstractMap) = MappedDomain(domain(d), map*mapping(d))

(*)(map::AbstractMap, domain::AbstractDomain) = apply_map(domain, map)

(*){N,T}(domain::AbstractDomain{N}, a::T) = scaling_map(a*ones(SVector{N,T})...) * domain

(+){N,T}(d::AbstractDomain{N}, x::SVector{N,T}) = AffineMap(eye(SMatrix{N,N,T}),x) * d
(+){N}(d::AbstractDomain{N}, x::AbstractVector) = d + SVector{N}(x)

show(io::IO, d::MappedDomain) =  print(io, "A mapped domain based on ", domain(d))
