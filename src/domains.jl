# domains.jl

abstract AbstractDomain{N,T}

dim{N,T}(::Type{AbstractDomain{N,T}}) = N
dim{D <: AbstractDomain}(::Type{D}) = dim(super(D))
dim(d::AbstractDomain) = dim(typeof(d))


numtype{N,T}(::Type{AbstractDomain{N,T}}) = T
numtype{D <: AbstractDomain}(::Type{D}) = numtype(super(D))
numtype(d::AbstractDomain) = numtype(typeof(d))

typealias AbstractDomain1d{T <: AbstractFloat} AbstractDomain{1,T}
typealias AbstractDomain2d{T <: AbstractFloat} AbstractDomain{2,T}
typealias AbstractDomain3d{T <: AbstractFloat} AbstractDomain{3,T}
typealias AbstractDomain4d{T <: AbstractFloat} AbstractDomain{4,T}

# Domains are evaluated using vectors to specify the points, except in 1D
# Provide fallback routine for users not using vectors in 1d
in{T,S <: Number}(x::S, d::AbstractDomain1d{T}) = in(Vec{1,S}(x), d)

# Check whether a value is in an interval, up to 10 times machine precision
in{T <: AbstractFloat, S <: Number}(x::S, a::T, b::T) = (a-10eps(T) <= x <= b+10eps(T))

# Fallback routine when evaluated on a grid. This routine is general, in order to avoid ambiguity
# with other routines later on. Dispatch on dimension is done by a different routine evalgrid below.
in(g::AbstractGrid, d::AbstractDomain) = evalgrid(g, d)

# left and right of domains falls back to bounding box domains
left(d::AbstractDomain) = left(boundingbox(d))
right(d::AbstractDomain) = right(boundingbox(d))

left(d::AbstractDomain, index::Int) = left(boundingbox(d), index)
right(d::AbstractDomain, index::Int) = right(boundingbox(d), index)

#in{N}(m::NTuple{N}, d::AbstractDomain) = in(Grid(boundingbox(d), m), d)

# Default methods for evaluation on a grid: the default is to call eval on the domain with 
# points as arguments. Domains that have faster grid evaluation routines may define their own version.
function evalgrid{N}(g::AbstractGrid{N}, d::AbstractDomain{N})
    z = zeros(Bool, size(g))
    for i in eachindex(g)
        z[i] = in(g[i], d)
    end
    z
end

show{N}(io::IO, v::Vec{N}) = print(io, Vector(v))

###############################################################################################
### An empty domain
###############################################################################################

immutable EmptyDomain{N,T} <: AbstractDomain{N,T}
end

EmptyDomain{T}(n::Int, ::Type{T}) = EmptyDomain{n,T}()


in(x::AnyVector, d::EmptyDomain) = false

# Arithmetic operations

(+)(d::EmptyDomain, x::Number) = d
(+)(x::Number, d::EmptyDomain) = d

(*)(d::EmptyDomain, x::Number) = d
(*)(x::Number, d::EmptyDomain) = d

(/)(d::EmptyDomain, x::Number) = d

(==)(d1::EmptyDomain, d2::EmptyDomain) = true

show(io::IO, d::EmptyDomain) = print(io, "an empty domain")


###############################################################################################
### The space R^n
###############################################################################################

immutable RnDomain{N,T} <: AbstractDomain{N,T}
end

RnDomain{T}(n::Int, ::Type{T}) = RnDomain{n,T}()

in(x::AnyVector, d::RnDomain) = true

# Arithmetic operations

(+)(d::RnDomain, x::Number) = d
(+)(x::Number, d::RnDomain) = d

(*)(d::RnDomain, x::Number) = d
(*)(x::Number, d::RnDomain) = d

(/)(d::RnDomain, x::Number) = d

(==)(d1::RnDomain, d2::RnDomain) = true

show(io::IO, e::RnDomain) = print(io, "the ", N, "-dimensional Euclidean space")


###############################################################################################
### An interval
###############################################################################################

immutable Interval{T} <: AbstractDomain{1,T}
    a     ::  T
    b     ::  T

    Interval(a = -one(T), b = one(T)) = new(a,b)
end

Interval() = Interval{Float64}()

Interval{T}(::Type{T}) = Interval{T}()

Interval{T}(a::T, b::T) = Interval{T}(a, b)

Interval{T <: AbstractFloat}(::Type{T}) = Interval(-one(T), one(T))

in(x::AnyVector, d::Interval) = in(x[1], d.a, d.b)

left(d::Interval) = d.a
right(d::Interval) = d.b

# Arithmetic operations

(+)(d::Interval, x::Number) = Interval(d.a+x, d.b+x)
(+)(x::Number, d::Interval) = d+x

(*)(d::Interval, x::Number) = Interval(x*d.a, x*d.b)
(*)(x::Number, d::Interval) = d*x

(/)(d::Interval, x::Number) = d * (1/x)

(==)(d1::Interval,d2::Interval) = (d1.a == d2.a) && (d1.b == d2.b)

boundingbox(d::Interval) = BBox(left(d), right(d))

show(io::IO, d::Interval) = print(io, "the interval [", d.a, ", ", d.b, "]")

const unitinterval = Interval()


###############################################################################################
### A circle
###############################################################################################

immutable Disk{T} <: AbstractDomain{2,T}
    radius    ::  T
    center    ::  Vec{2,T}

    Disk(radius = one(T), center = Vec{2,T}(0, 0)) = new(radius, center)
end

Disk() = Disk{Float64}()
Disk{T}(::Type{T}) = Disk{T}()
Disk{T}(radius::T) = Disk{T}(radius)
Disk{T}(radius::T, center::Vec{2,T}) = Disk{T}(radius, center)
Disk{T}(radius::T, center::Vector{T}) = Disk{T}(radius, center)

in{T}(x::AnyVector, c::Disk{T}) = (x[1]-c.center[1])^2 + (x[2]-c.center[2])^2 <= c.radius^2+10eps(T)

## Arithmetic operations

(+)(c::Disk, x::AnyVector) = Disk(c.radius, c.center+x)
(+)(x::AnyVector, c::Disk) = c+x

(*)(c::Disk, x::Number) = Disk(c.radius*x, c.center*x)
(*)(x::Number, c::Disk) = c*x

(/)(c::Disk, x::Number) = c * (1/x)

(==)(c1::Disk,c2::Disk) = (c1.radius == c2.radius) && (c1.center == c2.center)

boundingbox(c::Disk) = BBox((c.center[1]-c.radius,c.center[2]-c.radius),(c.center[1]+c.radius,c.center[2]+c.radius))

show(io::IO, c::Disk) = print(io, "a circle of radius ", c.radius, " centered at ", c.center)

const unitdisk = Disk()



###############################################################################################
### A 3D ball
###############################################################################################

immutable Ball{T} <: AbstractDomain{3,T}
    radius    ::  T
    center    ::  Vec{3,T}

    Ball(radius = one(T), center = Vec{3,T}(0, 0, 0)) = new(radius, center)
end

Ball() = Ball{Float64}()
Ball{T}(::Type{T}) = Ball{T}()
Ball{T}(radius::T) = Ball{T}(radius)
Ball{T}(radius::T, center::Vec{3,T}) = Ball{T}(radius, center)
Ball{T}(radius::T, center::Vector{T}) = Ball{T}(radius, center)

in(x::AnyVector, s::Ball) = (x[1]-s.center[1])^2 + (x[2]-s.center[2])^2 + (x[3]-s.center[3])^2 <= s.radius^2

## Arithmetic operations

(+)(s::Ball, x::AnyVector) = Ball(s.radius, s.center+x)
(+)(x::AnyVector, s::Ball) = Ball(s.radius, s.center+x)

(*)(s::Ball, x::Number) = Ball(s.radius * x, s.center * x)
(*)(x::Number, s::Ball) = s*x

(/)(s::Ball, x::Number) = s * (1/x)

(==)(s1::Ball, s2::Ball) = (s1.radius == s2.radius) && (s1.center == s2.center)

boundingbox(c::Ball) = BBox((c.center[1]-c.radius,c.center[2]-c.radius,c.center[3]-c.radius),(c.center[1]+c.radius,c.center[2]+c.radius,c.center[3]+c.radius))

show(io::IO, s::Ball) = print(io, "a sphere of radius ", s.radius, " centered at ", s.center)

const unitsphere = Ball()



###############################################################################################
### A Tensor Product of Domains
###############################################################################################

"""
A TensorProductDomain represents the tensor product of other domains.

immutable TensorProductDomain{TD,DN,LEN,N,T} <: AbstractDomain{N,T}

Parameters:
- TD is a tuple of (domain) types.
- DN is a tuple of the dimensions of each of the domains.
- LEN is the length of TG and GN
- N and T are the total dimension and numeric type of this grid.
"""
immutable TensorProductDomain{TD,DN,LEN,N,T} <: AbstractDomain{N,T}
	domains	::	TD

	TensorProductDomain(domains::Tuple) = new(domains)
end

tp_length{TD,DN,LEN,N,T}(::Type{TensorProductDomain{TD,DN,LEN,N,T}}) = LEN
tp_length(d::TensorProductDomain) = tp_length(typeof(d))

function TensorProductDomain(domains...)
    TD = typeof(domains)
    DN = map(dim, domains)
    LEN = length(domains)
    N = sum(DN)
    T = numtype(domains[1])
    TensorProductDomain{TD,DN,LEN,N,T}(domains)
end

⊗(d1::AbstractDomain, d2::AbstractDomain) = TensorProductDomain(d1, d2)
⊗(d1::TensorProductDomain, d2::TensorProductDomain) = TensorProductDomain(domainlist(d1)..., domainlist(d2)...)
⊗(d1::TensorProductDomain, d2::AbstractDomain) = TensorProductDomain(domainlist(d1)..., d2)
⊗(d1::AbstractDomain, d2::TensorProductDomain) = TensorProductDomain(d1, domainlist(d2)...)

tensorproduct(d::AbstractDomain, n) = TensorProductDomain([d for i=1:n]...)

subdomain(t::TensorProductDomain, i::Int) = t.domains[i]
domainlist(t::TensorProductDomain) = t.domains

# TODO: make this code for in more general!
# The problem is you can't slice a Vec, so the implementation below for AbstractArray does not work
# for FixedSizeArray's.
# All implementations below allocate memory (arrays) except the first one.
function in{TD,DN,N,T}(x::Vec{N,T}, t::TensorProductDomain{TD,DN,N,N,T})
    z1 = true
    for i = 1:N
        z1 = z1 & in(x[i], t.domains[i])
    end
    z1
end

function in{TD,DN,T}(x::Vec{3,T}, t::TensorProductDomain{TD,DN,2,3,T})
    N1 = DN[1]
    N2 = DN[2]
    d1 = subdomain(t, 1)
    d2 = subdomain(t, 2)
    x1 = Vec{N1,T}([x[j] for j=1:N1])
    x2 = Vec{N2,T}([x[j] for j=N1+1:N1+N2])
    in(x1, d1) && in(x2, d2)
end

function in{TD,DN,T}(x::Vec{4,T}, t::TensorProductDomain{TD,DN,2,4,T})
    N1 = DN[1]
    N2 = DN[2]
    d1 = subdomain(t, 1)
    d2 = subdomain(t, 2)
    x1 = Vec{N1,T}([x[j] for j=1:N1])
    x2 = Vec{N2,T}([x[j] for j=N1+1:N1+N2])
    in(x1, d1) && in(x2, d2)
end

function in{TD,DN,T}(x::Vec{4,T}, t::TensorProductDomain{TD,DN,3,4,T})
    d1 = subdomain(t, 1)
    d2 = subdomain(t, 2)
    d3 = subdomain(t, 3)
    x1 = Vec{N1,T}([x[j] for j=1:N1])
    x2 = Vec{N2,T}([x[j] for j=N1+1:N1+N2])
    x3 = Vec{N3,T}([x[j] for j=N1+N2+1:N1+N2+N3])
    in(x1, d1) && in(x2, d2) && in(x3, d3)
end


function in{TD,DN,LEN,N,T}(x::AbstractArray, t::TensorProductDomain{TD,DN,LEN,N,T})
    dc = 1
    z1 = true
    for i = 1:LEN
        z2 = in(x[dc:dc+DN[i]-1],t.domains[i])
        z1 = z1 & z2
        dc+=DN[i]
    end
    z1
end

# TODO: provide implementation of in for tensorproductgrids
 

(==)(t1::TensorProductDomain, t2::TensorProductDomain) = t1.domains==t2.domains

function boundingbox{TD,DN,LEN,N,T}(t::TensorProductDomain{TD,DN,LEN,N,T})
    dc = 1
    l = zeros(N)
    r = zeros(N)
    for i = 1:LEN
        for j = 1:DN[i]
            l[dc+j-1] = left(boundingbox(t.domains[i]), j)
            r[dc+j-1] = right(boundingbox(t.domains[i]), j)
        end
        dc += DN[i]
    end
    return BBox{N,T}(l, r)
end
 

function show(io::IO, t::TensorProductDomain)
    L = tp_length(t)
    for i = 1:L-1
        show(io, domainlist(t)[i])
        print(io, " x ")
    end
    show(io, domainlist(t)[L])
end



###############################################################################################
### An n-dimensional cube
###############################################################################################


Cube() = Cube(Val{3})
Cube{N}(::Type{Val{N}}) = Cube(Val{N}, Float64)

Cube{T}(::Type{Val{1}}, ::Type{T}) = Interval{T}()
Cube{N,T}(::Type{Val{N}}, ::Type{T}) = Interval{T}() ⊗ Cube(Val{N-1}, T)

Cube{T}(left::Vec{1,T}, right::Vec{1,T}) = Interval(left[1], right[1])
Cube{T}(left::Vec{2,T}, right::Vec{2,T}) = Interval(left[1], right[1]) ⊗ Interval(left[2], right[2])
Cube{T}(left::Vec{3,T}, right::Vec{3,T}) =
  Interval(left[1], right[1]) ⊗ Interval(left[2], right[2]) ⊗ Interval(left[3], right[3])
Cube{T}(left::Vec{4,T}, right::Vec{4,T}) =
  Interval(left[1], right[1]) ⊗ Interval(left[2], right[2]) ⊗ Interval(left[3], right[3]) ⊗ Interval(left[4], right[4])

Cube(n::Int) = Cube(Val{n})
Cube{N,T}(left::NTuple{N,T}, right::NTuple{N,T}) = Cube(Vec{N,T}(left), Vec{N,T}(right))

rectangle(a, b, c, d) = Interval(a,b) ⊗ Interval(c,d)

cube(a, b, c, d, e, f) = Interval(a,b) ⊗ Interval(c,d) ⊗ Interval(d,e)


const unitsquare = Cube(Val{2})
const unitcube = Cube(Val{3})



###############################################################################################
### A cylinder
###############################################################################################


Cylinder{T}(radius::T = one(T), length::T = one(T)) = Disk(radius) ⊗ Interval(zero(T),length)



###############################################################################################
### The union of two domains
###############################################################################################

# Type parameters N T D1 D2: dimension, numeric type, type of domain 1, type of domain 2.
# TODO: Make this a union of a list of domains
immutable DomainUnion{N,T,D1,D2} <: AbstractDomain{N,T}
    d1    ::  D1
    d2    ::  D2

    DomainUnion(d1::AbstractDomain{N,T}, d2::AbstractDomain{N,T}) = new(d1, d2)
end

DomainUnion{N,T}(d1::AbstractDomain{N,T}, d2::AbstractDomain{N,T}) = DomainUnion{N,T,typeof(d1),typeof(d2)}(d1, d2)

union(d1::AbstractDomain, d2::AbstractDomain) = (d1 == d2 ? d1 : DomainUnion(d1,d2))

function union(d1::Interval, d2::Interval)
    a = left(d1)
    b = right(d1)
    c = left(d2)
    d = right(d2)

    if (b < c) || (a > d)
        DomainUnion(d1, d2)
    else
        Interval(min(a, c), max(b, d))
    end
end


# The union of two domains corresponds to a logical OR of their characteristic functions
in(x::AnyVector, d::DomainUnion) = in(x, d.d1) || in(x, d.d2)

function in(g::AbstractGrid, d::DomainUnion)
    z1 = in(g, d.d1)
    z2 = in(g, d.d2)
    z1 | z2
end

(+)(d1::AbstractDomain, d2::AbstractDomain) = union(d1,d2)
(|)(d1::AbstractDomain, d2::AbstractDomain) = union(d1,d2)

(==)(d1::DomainUnion, d2::DomainUnion) = (d1.d1 == d2.d1) && (d1.d2 == d2.d2)

boundingbox(d::DomainUnion) = boundingbox(d.d1) ∪ boundingbox(d.d2)

function show(io::IO, d::DomainUnion)
    print(io, "A union of two domains: \n")
    print(io, "First domain: ", d.d1, "\n")
    print(io, "Second domain: ", d.d2, "\n")
end


###############################################################################################
### The intersection of two domains
###############################################################################################

immutable DomainIntersection{N,T,D1,D2} <: AbstractDomain{N,T}
    d1    ::  D1
    d2    ::  D2

    DomainIntersection(d1::AbstractDomain{N,T}, d2::AbstractDomain{N,T}) = new(d1, d2)
end

DomainIntersection{N,T}(d1::AbstractDomain{N,T},d2::AbstractDomain{N,T}) = DomainIntersection{N,T,typeof(d1),typeof(d2)}(d1, d2)

# The intersection of two domains corresponds to a logical AND of their characteristic functions
in(x::AnyVector, d::DomainIntersection) = in(x, d.d1) && in(x, d.d2)

function in(g::AbstractGrid, d::DomainIntersection)
    z1 = in(g, d.d1)
    z2 = in(g, d.d2)
    z1 & z2
end

(&)(d1::AbstractDomain, d2::AbstractDomain) = intersect(d1,d2)

intersect(d1::AbstractDomain, d2::AbstractDomain) = (d1 == d2 ? d1 : DomainIntersection(d1,d2))

function intersect{T}(d1::Interval{T}, d2::Interval{T})
    a = left(d1)
    b = right(d1)
    c = left(d2)
    d = right(d2)

    if (b < c) || (a > d)
        EmptyDomain{dim(d1),T}()
    else
        Interval(max(a, c), min(b, d))
    end
end

intersect{TD1,TD2,DN,LEN}(d1::TensorProductDomain{TD1,DN,LEN}, d2::TensorProductDomain{TD2,DN,LEN}) =
    TensorProductDomain([intersect(subdomain(d1,i), subdomain(d2,i)) for i in 1:LEN]...)


(==)(d1::DomainIntersection, d2::DomainIntersection) = (d1.d1 == d2.d1) && (d1.d2 == d2.d2)

boundingbox(d::DomainIntersection) = boundingbox(d.d1) ∩ boundingbox(d.d2)

function show(io::IO, d::DomainIntersection)
    print(io, "the intersection of two domains: \n")
    print(io, "    First domain: ", d.d1, "\n")
    print(io, "    Second domain: ", d.d1, "\n")
end


###############################################################################################
### The difference between two domains
###############################################################################################

immutable DomainDifference{N,T,D1,D2} <: AbstractDomain{N,T}
    d1    ::  D1
    d2    ::  D2

    DomainDifference(d1::AbstractDomain{N,T}, d2::AbstractDomain{N,T}) = new(d1, d2)
end

DomainDifference{N,T}(d1::AbstractDomain{N,T}, d2::AbstractDomain{N,T}) = DomainDifference{N,T,typeof(d1),typeof(d2)}(d1,d2)

# The difference between two domains corresponds to a logical AND NOT of their characteristic functions
in(x::AnyVector, d::DomainDifference) = in(x, d.d1) && (~in(x, d.d2))

function in(g::AbstractGrid, d::DomainDifference)
    z1 = in(g, d.d1)
    z2 = in(g, d.d2)
    z1 & (~z2)
end

(-)(d1::AbstractDomain, d2::AbstractDomain) = DomainDifference(d1,d2)
(\)(d1::AbstractDomain, d2::AbstractDomain) = DomainDifference(d1,d2)

(==)(d1::DomainDifference, d2::DomainDifference) = (d1.d1 == d2.d1) && (d1.d2 == d2.d2)

boundingbox(d::DomainDifference) = boundingbox(d.d1)

function show(io::IO, d::DomainDifference)
    print(io, "the difference of two domains: \n")
    print(io, "    First domain: ", d.d1, "\n")
    print(io, "    Second domain: ", d.d1, "\n")
end


###############################################################################################
### A revolved domain is a 2D-domain rotated about the X-axis
###############################################################################################

immutable RevolvedDomain{T,D} <: AbstractDomain{3,T}
    d     ::  D
end

revolve{T}(d::AbstractDomain{2,T}) = RevolvedDomain{T,typeof(d)}(d)

function in(x::AnyVector, d::RevolvedDomain)
    r = sqrt(x[2]^2+x[3])
    phi = atan2(x[2]/x[1])
    theta = acos(x[3]/r)
    in((x[1],r), d.d)
end

(==)(d1::RevolvedDomain, d2::RevolvedDomain) = (d1.d == d2.d)

boundingbox(d::RevolvedDomain) = BBox((left(d.d)[1],left(d.d)...),(right(d.d)[1],right(d.d)...))

function show(io::IO, r::RevolvedDomain)
    print(io, "the revolution of: ", r.d1)
end


###############################################################################################
### A rotated domain
###############################################################################################

immutable RotatedDomain{N,T,D} <: AbstractDomain{N,T}
    d                 ::  D
    angle             ::  Vector{T}
    rotationmatrix    ::  Array{T,2}

    # RotatedDomain(d,angle,rotationmatrix,box) = new(d, angle, rotationmatrix, box)
end

# Rotation in positive direction
rotationmatrix(theta) = Matrix2x2([cos(theta) -sin(theta); sin(theta) cos(theta)])
# Rotation about X-axis (phi), Y-axis (theta) and Z-axis (psi)
rotationmatrix(phi,theta,psi) = [cos(theta)*cos(psi) cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi) sin(phi)*sin(psi)-cos(phi)*sin(theta)*cos(psi); -cos(theta)*sin(psi) cos(phi)*cos(psi)-sin(phi)*sin(theta)*sin(psi) sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi); sin(theta) -sin(phi)*cos(theta) cos(phi)*cos(theta)]

RotatedDomain{T}(d::AbstractDomain{2,T}, theta) = RotatedDomain{2,T,typeof(d)}(d, [theta], rotationmatrix(theta))
# types annotated to remove ambiguity
RotatedDomain{T,D}(d::D, phi::T, theta::T, psi::T) = RotatedDomain{3,T,D}(d, [phi,theta,psi], rotationmatrix(phi,theta,psi))

rotate{T}(d::AbstractDomain{2,T}, theta) = RotatedDomain{2,T,typeof(d)}(d, theta)
rotate{T}(d::AbstractDomain{3,T}, phi::T, theta::T, psi::T) = RotatedDomain(d, phi, theta, psi)

in(x::AnyVector, d::RotatedDomain) = in(d.rotationmatrix*x, d.d)

(==)(d1::RotatedDomain, d2::RotatedDomain) = (d1.d == d2.d) && (d1.angle == d2.angle) #&& (d1.rotationmatrix == d2.rotationmatrix)

# very crude bounding box (doesn't work!!!)
boundingbox(r::RotatedDomain)= boundingbox(r.d)



###############################################################################################
### A scaled domain
###############################################################################################

# Note that the Euclidean plane is scaled, not just the domain itself.
# So the location of the origin matters. Two times a circle of radius 1 at a distance d of the origin
# becomes a circle of radius 2 at a distance 2d of the origin.
immutable ScaledDomain{N,T,D} <: AbstractDomain{N,T}
    domain      ::  D
    scalefactor ::  T

    ScaledDomain(domain::AbstractDomain{N,T}, scalefactor) = new(domain, scalefactor)
end

ScaledDomain{N,T}(domain::AbstractDomain{N,T}, scalefactor) = ScaledDomain{N,T,typeof(domain)}(domain, scalefactor)

domain(s::ScaledDomain) = d.domain

scalefactor(s::ScaledDomain) = d.scalefactor

function in(x::AnyVector, d::ScaledDomain)
    in(x/d.scalefactor, d.domain)
end

(*)(a::Number, d::AbstractDomain) = ScaledDomain(d, a)
(*)(a::Number, d::ScaledDomain) = ScaledDomain(domain(d), a*scalefactor(d))
(*)(d::AbstractDomain, a::Number) = a*d

boundingbox(s::ScaledDomain) = s.scalefactor * boundingbox(s.domain)



###############################################################################################
### A translated domain
###############################################################################################

immutable TranslatedDomain{N,T,D} <: AbstractDomain{N,T}
    domain  ::  D
    trans   ::  Vec{N,T}

    TranslatedDomain(domain::AbstractDomain{N,T}, trans) = new(domain, trans)
end

TranslatedDomain{N,T}(domain::AbstractDomain{N,T}, trans::AnyVector) = TranslatedDomain{N,T,typeof(domain)}(domain, trans)

domain(d::TranslatedDomain) = d.domain

translationvector(d::TranslatedDomain) = d.trans

function in(x::AnyVector, d::TranslatedDomain)
    in(x-d.trans, d.domain)
end

(+)(d::AbstractDomain, trans::AnyVector) = TranslatedDomain(d, trans)
(+)(d::TranslatedDomain, trans::AnyVector) = TranslatedDomain(domain(d), trans+translationvector(d))
(+)(trans::AnyVector, d::AbstractDomain) = d + a

boundingbox(d::TranslatedDomain) = boundingbox(domain(d)) + trans



###############################################################################################
### A collection of domains
###############################################################################################

type DomainCollection{N,T} <: AbstractDomain{N,T}
    list    ::  Vector{AbstractDomain{N,T}}
end

DomainCollection{N,T}(d::AbstractDomain{N,T}) = DomainCollection{N,T}([d])

length(d::DomainCollection) = length(d.list)

domain(d::DomainCollection, i) = d.list[i]

function in(x::AnyVector, d::DomainCollection)
    reduce( |, map( u -> in(x, u), d.list))
end

function in(g::AbstractGrid, d::DomainCollection)
    z1 = in(g, d.list[1])
    for i = 2:length(d.list)
        z2 = in(g, d.list[i])
        z1 = z1 | z2
    end
    z1
end

push!(dc::DomainCollection, d::AbstractDomain) = push!(dc.list, d)

(==)(d1::DomainCollection, d2::DomainCollection) = reduce(&, map( (x,y) -> x==y, d1.list, d2.list))

 function boundingbox(d::DomainCollection)
     ubox=boundingbox(d.list[1])
     for i = 2:length(d.list)
         ubox = union(ubox, boundingbox(d.list[1]))
     end
     ubox
 end

 
show(io::IO, d::DomainCollection) = print(io, "a collection of ", length(d.list), " domains")


 
##########################################################################
### Assorted Domains
##########################################################################


function randomcircles(n)
    list = [Disk(0.2, (2*rand(2)-1)*0.8) for i=1:n]
    DC = DomainCollection(list[1])
    for i = 2:n
        push!(DC.list, list[i])
    end
    DC.box = BBox(-1.0, 1.0, -1.0, 1.0)
    DC
end


###
# The atomium: a famous building in Belgium
###

function atomium()
    sphere1 = Ball(0.25)
    spheres = DomainCollection(sphere1)
    push!(spheres, sphere1 + [ 0.6, 0.6, 0.6])
    push!(spheres, sphere1 + [ 0.6, 0.6,-0.6])
    push!(spheres, sphere1 + [ 0.6,-0.6, 0.6])
    push!(spheres, sphere1 + [ 0.6,-0.6,-0.6])
    push!(spheres, sphere1 + [-0.6, 0.6, 0.6])
    push!(spheres, sphere1 + [-0.6, 0.6,-0.6])
    push!(spheres, sphere1 + [-0.6,-0.6, 0.6])
    push!(spheres, sphere1 + [-0.6,-0.6,-0.6])
    cyl1 = Cylinder(0.10, 1.2)
    push!(spheres, cyl1 + [-0.6, 0.6, 0.6]);
    push!(spheres, cyl1 + [-0.6,-0.6, 0.6]);
    push!(spheres, cyl1 + [-0.6, 0.6,-0.6]);
    push!(spheres, cyl1 + [-0.6,-0.6,-0.6]);
    cyl2 = rotate(cyl1, 0.0, 0.0, pi/2.0)
    push!(spheres, cyl2 + [ 0.6, -0.6, 0.6])
    push!(spheres, cyl2 + [-0.6, -0.6, 0.6])
    push!(spheres, cyl2 + [ 0.6, -0.6,-0.6])
    push!(spheres, cyl2 + [-0.6, -0.6,-0.6])
    cyl2b = rotate(cyl1, 0.0, pi/2.0, 0.0)
    push!(spheres, cyl2b + [ 0.6,  0.6, 0.6])
    push!(spheres, cyl2b + [-0.6,  0.6, 0.6])
    push!(spheres, cyl2b + [ 0.6, -0.6, 0.6])
    push!(spheres, cyl2b + [-0.6, -0.6, 0.6])
    cyl3 = Cylinder(0.10, 1.2*sqrt(3))
    cyl3 = rotate(cyl3, 0.0, asin(1/sqrt(3)), 0.0)
    cyl3 = rotate(cyl3, 0.0, 0.0, pi/4)
    push!(spheres, cyl3 + [ -0.6, -0.6, +0.6])
    cyl4 = Cylinder(0.10, 1.2*sqrt(3))
    cyl4 = rotate(cyl4, 0.0, -asin(1/sqrt(3)), 0.0)
    cyl4 = rotate(cyl4, 0.0, 0.0, pi/4)
    push!(spheres, cyl4 + [ -0.6, -0.6, -0.6])
    cyl5 = Cylinder(0.10, 1.2*sqrt(3))
    cyl5 = rotate(cyl5, 0.0, asin(1/sqrt(3)), 0.0)
    cyl5 = rotate(cyl5, 0.0, 0.0, -pi/4)
    push!(spheres, cyl5 + [ -0.6, +0.6, +0.6])
    cyl6 = Cylinder(0.10, 1.2*sqrt(3))
    cyl6 = rotate(cyl6, 0.0, -asin(1/sqrt(3)), 0.0)
    cyl6 = rotate(cyl6, 0.0, 0.0, -pi/4)
    push!(spheres, cyl6 + [ -0.6, +0.6, -0.6])
    spheres.box = unitbox3
    atomium = spheres
end


## An affinely mapped domain

# Find the map that maps a and b to c and d
# function affinemap(a, b, c, d)

# end

# type AffineMapDomain{N,T} <: AbstractDomain{N,T}
#     A
#     box
# end









