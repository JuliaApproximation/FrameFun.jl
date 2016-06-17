# derived_domains.jl


###############################################################################################
### A domain described by a characteristic function
###############################################################################################

immutable Characteristic{N,T} <: AbstractDomain{N}
    char    ::  Function
    box    ::  BBox{N,T}
end


in(x::Vec, c::Characteristic) = c.char(x)

boundingbox(c::Characteristic) = c.box

show(io::IO, c::Characteristic) = print(io, "a domain described by a characteristic function")




###############################################################################################
### The union of two domains
###############################################################################################

immutable DomainUnion{D1,D2,N} <: AbstractDomain{N}
    d1    ::  D1
    d2    ::  D2

    DomainUnion(d1::AbstractDomain{N}, d2::AbstractDomain{N}) = new(d1, d2)
end

DomainUnion{N}(d1::AbstractDomain{N}, d2::AbstractDomain{N}) = DomainUnion{typeof(d1),typeof(d2),N}(d1, d2)

union(d1::AbstractDomain, d2::AbstractDomain) = (d1 == d2 ? d1 : DomainUnion(d1,d2))


# The union of two domains corresponds to a logical OR of their characteristic functions
in(x::Vec, d::DomainUnion) = in(x, d.d1) || in(x, d.d2)

function in(g::AbstractGrid, d::DomainUnion)
    z1 = in(g, d.d1)
    z2 = in(g, d.d2)
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


###############################################################################################
### The intersection of two domains
###############################################################################################

immutable DomainIntersection{D1,D2,N} <: AbstractDomain{N}
    d1    ::  D1
    d2    ::  D2

    DomainIntersection(d1::AbstractDomain{N}, d2::AbstractDomain{N}) = new(d1, d2)
end

DomainIntersection{N}(d1::AbstractDomain{N},d2::AbstractDomain{N}) = DomainIntersection{typeof(d1),typeof(d2),N}(d1, d2)

# The intersection of two domains corresponds to a logical AND of their characteristic functions
in(x::Vec, d::DomainIntersection) = in(x, d.d1) && in(x, d.d2)

function in(g::AbstractGrid, d::DomainIntersection)
    z1 = in(g, d.d1)
    z2 = in(g, d.d2)
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


###############################################################################################
### The difference between two domains
###############################################################################################

immutable DomainDifference{D1,D2,N} <: AbstractDomain{N}
    d1    ::  D1
    d2    ::  D2

    DomainDifference(d1::AbstractDomain{N}, d2::AbstractDomain{N}) = new(d1, d2)
end

DomainDifference{N}(d1::AbstractDomain{N}, d2::AbstractDomain{N}) = DomainDifference{typeof(d1),typeof(d2),N}(d1,d2)

setdiff(d1::AbstractDomain, d2::AbstractDomain) = DomainDifference(d1, d2)

# The difference between two domains corresponds to a logical AND NOT of their characteristic functions
in(x::Vec, d::DomainDifference) = in(x, d.d1) && (~in(x, d.d2))

function in(g::AbstractGrid, d::DomainDifference)
    z1 = in(g, d.d1)
    z2 = in(g, d.d2)
    z1 & (~z2)
end

(-)(d1::AbstractDomain, d2::AbstractDomain) = setdiff(d1, d2)
(\)(d1::AbstractDomain, d2::AbstractDomain) = setdiff(d1, d2)


boundingbox(d::DomainDifference) = boundingbox(d.d1)

function show(io::IO, d::DomainDifference)
    print(io, "the difference of two domains: \n")
    print(io, "    First domain: ", d.d1, "\n")
    print(io, "    Second domain: ", d.d2, "\n")
end


###############################################################################################
### A revolved domain is a 2D-domain rotated about the X-axis
###############################################################################################

immutable RevolvedDomain{D} <: AbstractDomain{3}
    d     ::  D
end

revolve(d::AbstractDomain{2}) = RevolvedDomain(d)

function in(x::Vec, d::RevolvedDomain)
    r = sqrt(x[2]^2+x[3])
    phi = atan2(x[2]/x[1])
    theta = acos(x[3]/r)
    in((x[1],r), d.d)
end


boundingbox(d::RevolvedDomain) = BBox((left(d.d)[1],left(d.d)...),(right(d.d)[1],right(d.d)...))

function show(io::IO, r::RevolvedDomain)
    print(io, "the revolution of: ", r.d1)
end


###############################################################################################
### A rotated domain
###############################################################################################

immutable RotatedDomain{D,T,N} <: AbstractDomain{N}
   d                 ::  D
   angle             ::  Vector{T}
   rotationmatrix    ::  Mat{N,N,T}

    RotatedDomain(d,angle,rotationmatrix) = new(d, angle, rotationmatrix)
end

# Rotation in positive direction
rotationmatrix(theta) = Mat{2,2,typeof(theta)}([cos(theta) -sin(theta); sin(theta) cos(theta)])
# Rotation about X-axis (phi), Y-axis (theta) and Z-axis (psi)
rotationmatrix(phi,theta,psi) =
   Mat{3,3,typeof(phi)}([cos(theta)*cos(psi) cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi) sin(phi)*sin(psi)-cos(phi)*sin(theta)*cos(psi); -cos(theta)*sin(psi) cos(phi)*cos(psi)-sin(phi)*sin(theta)*sin(psi) sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi); sin(theta) -sin(phi)*cos(theta) cos(phi)*cos(theta)])

RotatedDomain{N,T}(d::AbstractDomain{N}, angle::Vector{T}, m::Mat{N,N,T} = rotationmatrix(theta)) =
   RotatedDomain{typeof(d),T,N}(d, angle, m)

RotatedDomain(d::AbstractDomain{2}, theta::Number) = RotatedDomain{2,typeof(theta),typeof(d)}(d, [theta], rotationmatrix(theta))
# types annotated to remove ambiguity
RotatedDomain{T,D}(d::D, phi::T, theta::T, psi::T) = RotatedDomain{3,T,D}(d, [phi,theta,psi], rotationmatrix(phi,theta,psi))

rotate{T}(d::AbstractDomain{2}, theta::T) = RotatedDomain(d, theta)
rotate{T}(d::AbstractDomain{3}, phi::T, theta::T, psi::T) = RotatedDomain(d, phi, theta, psi)

in(x::Vec, d::RotatedDomain) = in(d.rotationmatrix*x, d.d
)
(==)(d1::RotatedDomain, d2::RotatedDomain) = (d1.d == d2.d) && (d1.angle == d2.angle) #&& (d1.rotationmatrix == d2.rotationmatrix)

# very crude bounding box (doesn't work!!!)
boundingbox(r::RotatedDomain)= sqrt(2)*boundingbox(r.d)



###############################################################################################
### A scaled domain
###############################################################################################

# Note that the Euclidean plane is scaled, not just the domain itself.
# Hence, the location of the origin matters. Two times a circle of radius 1 at a distance d of the origin
# becomes a circle of radius 2 at a distance 2d of the origin.
immutable ScaledDomain{D,T,N} <: AbstractDomain{N}
    domain      ::  D
    scalefactor ::  T

    ScaledDomain(domain::AbstractDomain{N}, scalefactor) = new(domain, scalefactor)
end

ScaledDomain{N,T}(domain::AbstractDomain{N}, scalefactor::T) = ScaledDomain{typeof(domain),T,N}(domain, scalefactor)

domain(s::ScaledDomain) = d.domain

scalefactor(s::ScaledDomain) = d.scalefactor

function in(x::Vec, d::ScaledDomain)
    in(x/d.scalefactor, d.domain)
end

(*)(d::AbstractDomain, x::Number) = ScaledDomain(d, x)
(*)(d::ScaledDomain, x::Number) = ScaledDomain(domain(d), x*scalefactor(d))

boundingbox(s::ScaledDomain) = scalefactor(s) * boundingbox(s.domain)



###############################################################################################
### A translated domain
###############################################################################################

immutable TranslatedDomain{D,T,N} <: AbstractDomain{N}
    domain  ::  D
    trans   ::  Vec{N,T}

    TranslatedDomain(domain::AbstractDomain{N}, trans) = new(domain, trans)
end

TranslatedDomain{N}(domain::AbstractDomain{N}, trans::Vec{N}) = TranslatedDomain{typeof(domain),eltype(trans),N}(domain, trans)

domain(d::TranslatedDomain) = d.domain

translationvector(d::TranslatedDomain) = d.trans

function in(x::Vec, d::TranslatedDomain)
    in(x-d.trans, d.domain)
end

(+)(d::AbstractDomain, trans::Vec) = TranslatedDomain(d, trans)
(+)(d::TranslatedDomain, trans::Vec) = TranslatedDomain(domain(d), trans+translationvector(d))

boundingbox(d::TranslatedDomain) = boundingbox(domain(d)) + translationvector(d)



###############################################################################################
### A collection of domains
###############################################################################################

type DomainCollection{N} <: AbstractDomain{N}
    list    ::  Vector{AbstractDomain{N}}
end

DomainCollection(d::AbstractDomain) = DomainCollection([d])


length(d::DomainCollection) = length(d.list)

domain(d::DomainCollection, i) = d.list[i]

# Iteration over the domain list
start(d::DomainCollection) = start(d.list)

next(d::DomainCollection, state) = next(d.list, state)

done(d::DomainCollection, state) = done(d.list, state)

function in(x::Vec, dc::DomainCollection)
    z = false
    for d in dc
        z = z || in(x, d)
    end
    z
end

function evalgrid!(z, g::AbstractGrid, dc::DomainCollection)
    for d in dc
        evalgrid!(z, g, d)
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



## An affinely mapped domain

# Find the map that maps a and b to c and d
# function affinemap(a, b, c, d)

# end

# type AffineMapDomain{N,T} <: AbstractDomain{N,T}
#     A
#     box
# end
