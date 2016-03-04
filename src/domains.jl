# domains.jl

"An N-dimensional domain."
abstract AbstractDomain{N}

dim{N}(d::AbstractDomain{N}) = N


typealias AbstractDomain1d AbstractDomain{1}
typealias AbstractDomain2d AbstractDomain{2}
typealias AbstractDomain3d AbstractDomain{3}
typealias AbstractDomain4d AbstractDomain{4}

# left and right of domains falls back to bounding box domains
left(d::AbstractDomain) = left(boundingbox(d))
right(d::AbstractDomain) = right(boundingbox(d))

left(d::AbstractDomain, i::Int) = left(boundingbox(d), i)
right(d::AbstractDomain, i::Int) = right(boundingbox(d), i)

# Domains are evaluated using vectors to specify the points, except in 1D
# Provide fallback routine for users not using vectors in 1d
in(x::Number, d::AbstractDomain1d) = in(Vec(x), d)

# Check whether a value is in an interval, up to 10 times machine precision
in{T <: AbstractFloat}(x::Number, a::T, b::T) = (a-10eps(T) <= x <= b+10eps(T))
in{T <: Number}(x::Number, a::T, b::T) = a <= x <= b

# Evaluation on a grid should be implemented by evalgrid for each domain
in{N}(g::AbstractGrid{N}, d::AbstractDomain{N}) = evalgrid(g, d)


# Default methods for evaluation on a grid: the default is to call eval on the domain with 
# points as arguments. Domains that have faster grid evaluation routines may define their own version.
evalgrid(g::AbstractGrid, d::AbstractDomain) = evalgrid!(zeros(Bool, size(g)), g, d)

# Note that evalgrid! only updates the result - it should be initialized to all false!
# The idea is that you can chain different calls to evalgrid (as used in DomainCollection below)
function evalgrid!{N}(result, g::AbstractGrid{N}, d::AbstractDomain{N})
    for i in eachindex(g)
        result[i] |= in(g[i], d)
    end
    result
end

show{N}(io::IO, v::Vec{N}) = print(io, Vector(v))


## Arithmetics

# Make suredomains only need to implement addition/multiplication with numbers to the right
(+)(x::Number, d::AbstractDomain) = d + x
(+)(x::AnyVector, d::AbstractDomain) = d + x

(*)(x::Number, d::AbstractDomain) = d * x

(/)(d::AbstractDomain, x::Number) = d * (1/x)


###############################################################################################
### An empty domain
###############################################################################################

immutable EmptyDomain{N} <: AbstractDomain{N}
end

EmptyDomain{N}(::Type{Val{N}} = Val{1}) = EmptyDomain{N}()
EmptyDomain(n::Int) = EmptyDomain(Val{n})

in(x::AnyVector, d::EmptyDomain) = false

# Arithmetic operations

(+)(d::EmptyDomain, x::Number) = d

(*)(d::EmptyDomain, x::Number) = d


show(io::IO, d::EmptyDomain) = print(io, "an empty domain")


###############################################################################################
### The space R^n
###############################################################################################

immutable RnDomain{N} <: AbstractDomain{N}
end

RnDomain{N}(::Type{Val{N}} = Val{1}) = RnDomain{N}()
RnDomain(n::Int) = RnDomain(Val{n})

in(x::AnyVector, d::RnDomain) = true

# Arithmetic operations

(+)(d::RnDomain, x::Number) = d

(*)(d::RnDomain, x::Number) = d


show(io::IO, e::RnDomain) = print(io, "the ", N, "-dimensional Euclidean space")


###############################################################################################
### An interval
###############################################################################################

immutable Interval{T} <: AbstractDomain{1}
    a     ::  T
    b     ::  T
end

Interval() = Interval(-1, 1)

Interval{T}(::Type{T}) = Interval{T}(-1, 1)


in(x::AnyVector, d::Interval) = in(x[1], d.a, d.b)

left(d::Interval) = d.a
right(d::Interval) = d.b

# Arithmetic operations

(+)(d::Interval, x::Number) = Interval(d.a+x, d.b+x)

(*)(d::Interval, x::Number) = Interval(x*d.a, x*d.b)


boundingbox(d::Interval) = BBox(left(d), right(d))

show(io::IO, d::Interval) = print(io, "the interval [", d.a, ", ", d.b, "]")

const unitinterval = Interval()


###############################################################################################
### A circle
###############################################################################################

immutable Disk{S,T} <: AbstractDomain{2}
    radius    ::  S
    center    ::  Vec{2,T}

    Disk(radius = one(S), center = Vec(0,0)) = new(radius, center)
end

Disk() = Disk{Int,Float64}()
Disk{T}(::Type{T}) = Disk{T,T}()

Disk{T}(radius::T) = Disk{T,T}(radius)
Disk{S,T}(radius::S, center::Vec{2,T}) = Disk{S,T}(radius, center)
Disk(radius, center::AbstractVector) = Disk(radius, Vec(center...))

in(x::AnyVector, c::Disk) = (x[1]-c.center[1])^2 + (x[2]-c.center[2])^2 <= c.radius^2

## Arithmetic operations

(+)(c::Disk, x::AnyVector) = Disk(c.radius, c.center+x)

(*)(c::Disk, x::Number) = Disk(c.radius*x, c.center*x)


boundingbox(c::Disk) = BBox((c.center[1]-c.radius,c.center[2]-c.radius),(c.center[1]+c.radius,c.center[2]+c.radius))

show(io::IO, c::Disk) = print(io, "a circle of radius ", c.radius, " centered at ", c.center)

const unitdisk = Disk()


###############################################################################################
### A domain described by a characteristic function
###############################################################################################

immutable Characteristic{N,T} <: AbstractDomain{N}
    char    ::  Function
    box    ::  BBox{N,T}
end


in(x::AnyVector, c::Characteristic) = c.char(x)

boundingbox(c::Characteristic) = c.box

show(io::IO, c::Characteristic) = print(io, "a domain described by a characteristic function")




###############################################################################################
### A 3D ball
###############################################################################################

immutable Ball{S,T} <: AbstractDomain{3}
    radius    ::  S
    center    ::  Vec{3,T}

    Ball(radius = one(S), center = Vec{3,T}(0, 0, 0)) = new(radius, center)
end

Ball() = Ball{Int,Float64}()
Ball{T}(::Type{T}) = Ball{T,T}()

Ball{T}(radius::T) = Ball{T,T}(radius)
Ball{S,T}(radius::S, center::Vec{3,T}) = Ball{S,T}(radius, center)
Ball(radius, center::AbstractVector) = Ball(radius, Vec(center...))


in(x::AnyVector, s::Ball) = (x[1]-s.center[1])^2 + (x[2]-s.center[2])^2 + (x[3]-s.center[3])^2 <= s.radius^2

## Arithmetic operations

(+)(s::Ball, x::AnyVector) = Ball(s.radius, s.center+x)

(*)(s::Ball, x::Number) = Ball(s.radius * x, s.center * x)


boundingbox(c::Ball) = BBox((c.center[1]-c.radius,c.center[2]-c.radius,c.center[3]-c.radius),(c.center[1]+c.radius,c.center[2]+c.radius,c.center[3]+c.radius))

show(io::IO, s::Ball) = print(io, "a sphere of radius ", s.radius, " centered at ", s.center)

const unitball = Ball()



###############################################################################################
### A Tensor Product of Domains
###############################################################################################

"""
A TensorProductDomain represents the tensor product of other domains.

immutable TensorProductDomain{TD,DN,LEN,N,T} <: AbstractDomain{N,T}

Parameters:
- TD is a tuple of (domain) types
- DN is a tuple of the dimensions of each of the domains
- LEN is the length of TG and GN
- N is the total dimension of this domain
"""
immutable TensorProductDomain{TD,DN,LEN,N} <: AbstractDomain{N}
	domains	::	TD

	TensorProductDomain(domains::Tuple) = new(domains)
end

tp_length{TD,DN,LEN,N}(d::TensorProductDomain{TD,DN,LEN,N}) = LEN

function TensorProductDomain(domains...)
    TD = typeof(domains)
    DN = map(dim, domains)
    LEN = length(domains)
    N = sum(DN)
    TensorProductDomain{TD,DN,LEN,N}(domains)
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
function in{TD,DN,N}(x::Vec{N}, t::TensorProductDomain{TD,DN,N,N})
    z1 = true
    for i = 1:N
        z1 = z1 & in(x[i], t.domains[i])
    end
    z1
end

function in{TD,DN,T}(x::Vec{3,T}, t::TensorProductDomain{TD,DN,2,3})
    N1 = DN[1]
    N2 = DN[2]
    d1 = subdomain(t, 1)
    d2 = subdomain(t, 2)
    x1 = Vec{N1,T}([x[j] for j=1:N1])
    x2 = Vec{N2,T}([x[j] for j=N1+1:N1+N2])
    in(x1, d1) && in(x2, d2)
end

function in{TD,DN,T}(x::Vec{4,T}, t::TensorProductDomain{TD,DN,2,4})
    N1 = DN[1]
    N2 = DN[2]
    d1 = subdomain(t, 1)
    d2 = subdomain(t, 2)
    x1 = Vec{N1,T}([x[j] for j=1:N1])
    x2 = Vec{N2,T}([x[j] for j=N1+1:N1+N2])
    in(x1, d1) && in(x2, d2)
end

function in{TD,DN,T}(x::Vec{4,T}, t::TensorProductDomain{TD,DN,3,4})
    d1 = subdomain(t, 1)
    d2 = subdomain(t, 2)
    d3 = subdomain(t, 3)
    x1 = Vec{N1,T}([x[j] for j=1:N1])
    x2 = Vec{N2,T}([x[j] for j=N1+1:N1+N2])
    x3 = Vec{N3,T}([x[j] for j=N1+N2+1:N1+N2+N3])
    in(x1, d1) && in(x2, d2) && in(x3, d3)
end


function in{TD,DN,LEN,N}(x::AbstractArray, t::TensorProductDomain{TD,DN,LEN,N})
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
 

function boundingbox{TD,DN,LEN,N}(t::TensorProductDomain{TD,DN,LEN,N})
    box = boundingbox(t.domains[1])
    for i = 2:LEN
        box = box ⊗ boundingbox(t.domains[i])
    end
    box
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

Cube() = Cube(Val{1})

Cube(::Type{Val{1}}) = Interval()
Cube{N}(::Type{Val{N}}) = Interval() ⊗ Cube(Val{N-1})
Cube(n::Int) = Cube(Val{n})

Cube(left::Number, right::Number) = Interval(left, right)
Cube(left, right) = TensorProductDomain(map(Interval, left, right)...)

rectangle(a, b, c, d) = Interval(a,b) ⊗ Interval(c,d)

cube(a, b, c, d, e, f) = Interval(a,b) ⊗ Interval(c,d) ⊗ Interval(e,f)


const unitsquare = Cube(Val{2})
const unitcube = Cube(Val{3})



###############################################################################################
### A cylinder
###############################################################################################


Cylinder(radius = 1, length = 1) = Disk(radius) ⊗ Interval(0,length)



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


boundingbox(d::DomainUnion) = boundingbox(d.d1) ∪ boundingbox(d.d2)

function show(io::IO, d::DomainUnion)
    print(io, "A union of two domains: \n")
    print(io, "First domain: ", d.d1, "\n")
    print(io, "Second domain: ", d.d2, "\n")
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
in(x::AnyVector, d::DomainIntersection) = in(x, d.d1) && in(x, d.d2)

function in(g::AbstractGrid, d::DomainIntersection)
    z1 = in(g, d.d1)
    z2 = in(g, d.d2)
    z1 & z2
end

(&)(d1::AbstractDomain, d2::AbstractDomain) = intersect(d1,d2)

intersect(d1::AbstractDomain, d2::AbstractDomain) = (d1 == d2 ? d1 : DomainIntersection(d1,d2))

function intersect(d1::Interval, d2::Interval)
    a = left(d1)
    b = right(d1)
    c = left(d2)
    d = right(d2)

    if (b < c) || (a > d)
        EmptyDomain(Val{1})
    else
        Interval(max(a, c), min(b, d))
    end
end

intersect{TD1,TD2,DN,LEN}(d1::TensorProductDomain{TD1,DN,LEN}, d2::TensorProductDomain{TD2,DN,LEN}) =
    TensorProductDomain([intersect(subdomain(d1,i), subdomain(d2,i)) for i in 1:LEN]...)



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
in(x::AnyVector, d::DomainDifference) = in(x, d.d1) && (~in(x, d.d2))

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

function in(x::AnyVector, d::RevolvedDomain)
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

#immutable RotatedDomain{D,T,N} <: AbstractDomain{N}
#    d                 ::  D
#    angle             ::  Vector{T}
#    rotationmatrix    ::  Array{T,2}

    # RotatedDomain(d,angle,rotationmatrix,box) = new(d, angle, rotationmatrix, box)
#end

# Rotation in positive direction
#rotationmatrix(theta) = Matrix2x2([cos(theta) -sin(theta); sin(theta) cos(theta)])
# Rotation about X-axis (phi), Y-axis (theta) and Z-axis (psi)
#rotationmatrix(phi,theta,psi) =
#    [cos(theta)*cos(psi) cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi) sin(phi)*sin(psi)-cos(phi)*sin(theta)*cos(psi); -cos(theta)*sin(psi) cos(phi)*cos(psi)-sin(phi)*sin(theta)*sin(psi) sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi); sin(theta) -sin(phi)*cos(theta) cos(phi)*cos(theta)]

#RotatedDomain{N}(d::AbstractDomain{N}, angle::Vector{T}, m::Array{T,2} = rotationmatrix(theta)) =
#    RotatedDomain{typeof(d),T,N}(d, angle, m)

#RotatedDomain(d::AbstractDomain{2}, theta) = RotatedDomain{2,T,typeof(d)}(d, [theta], rotationmatrix(theta))
# types annotated to remove ambiguity
#RotatedDomain{T,D}(d::D, phi::T, theta::T, psi::T) = RotatedDomain{3,T,D}(d, [phi,theta,psi], rotationmatrix(phi,theta,psi))

#rotate{T}(d::AbstractDomain{2,T}, theta) = RotatedDomain{2,T,typeof(d)}(d, theta)
#rotate{T}(d::AbstractDomain{3,T}, phi::T, theta::T, psi::T) = RotatedDomain(d, phi, theta, psi)

#in(x::AnyVector, d::RotatedDomain) = in(d.rotationmatrix*x, d.d)

#(==)(d1::RotatedDomain, d2::RotatedDomain) = (d1.d == d2.d) && (d1.angle == d2.angle) #&& (d1.rotationmatrix == d2.rotationmatrix)

# very crude bounding box (doesn't work!!!)
#boundingbox(r::RotatedDomain)= boundingbox(r.d)



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

function in(x::AnyVector, d::ScaledDomain)
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

TranslatedDomain{N}(domain::AbstractDomain{N}, trans::AnyVector) = TranslatedDomain{typeof(domain),eltype(trans),N}(domain, trans)

domain(d::TranslatedDomain) = d.domain

translationvector(d::TranslatedDomain) = d.trans

function in(x::AnyVector, d::TranslatedDomain)
    in(x-d.trans, d.domain)
end

(+)(d::AbstractDomain, trans::AnyVector) = TranslatedDomain(d, trans)
(+)(d::TranslatedDomain, trans::AnyVector) = TranslatedDomain(domain(d), trans+translationvector(d))

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

function in(x::AnyVector, dc::DomainCollection)
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


 
##########################################################################
### Assorted Domains
##########################################################################


function randomcircles(n)
    list = AbstractDomain[Disk(0.2, (2*rand(2)-1)*0.8) for i=1:n]
    DomainCollection{2}(list)
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









