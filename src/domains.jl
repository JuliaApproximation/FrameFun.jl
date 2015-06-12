# domains.jl

abstract AbstractDomain{N,T <: FloatingPoint}

dim{N,T}(::AbstractDomain{N,T}) = N
dim{N,T}(::Type{AbstractDomain{N,T}}) = N
dim{D <: AbstractDomain}(::Type{D}) = dim(super(D))

numtype{N,T}(::AbstractDomain{N,T}) = T
numtype{N,T}(::Type{AbstractDomain{N,T}}) = T
numtype{D <: AbstractDomain}(::Type{D}) = numtype(super(D))

typealias AbstractDomain1d{T <: FloatingPoint} AbstractDomain{1,T}
typealias AbstractDomain2d{T <: FloatingPoint} AbstractDomain{2,T}
typealias AbstractDomain3d{T <: FloatingPoint} AbstractDomain{3,T}
typealias AbstractDomain4d{T <: FloatingPoint} AbstractDomain{4,T}

# Domains are evaluated using vectors to specify the points, except in 1D
# Provide fallback routine for users not using vectors in 1d
in{T,S <: Number}(x::S, d::AbstractDomain1d{T}) = in([x], d)

# Check whether a value is in an interval, up to 10 times machine precision
in{T <: FloatingPoint, S <: Number}(x::S, a::T, b::T) = (a-10eps(T) <= x <= b+10eps(T))

# Fallback routine when evaluated on a grid. This routine is general, in order to avoid ambiguity
# with other routines later on. Dispatch on dimension is done by a different routine evalgrid below.
in(g::AbstractGrid, d::AbstractDomain) = evalgrid(g, d)

#in{N}(m::NTuple{N}, d::AbstractDomain) = in(Grid(box(d), m), d)

# Default methods for evaluation on a grid: the default is to call eval on the domain with 
# points as arguments. Domains that have faster grid evaluation routines may define their own version.
function evalgrid{N}(g::AbstractGrid{N}, d::AbstractDomain{N})
    z = zeros(Bool, size(g))
    for i in eachindex(g)
      z[i] = in(g[i], d)
    end
    z
end


###############################################################################################
### An empty domain
###############################################################################################

immutable EmptyDomain{N,T} <: AbstractDomain{N,T}
end

EmptyDomain{T}(n::Int, ::Type{T}) = EmptyDomain{n,T}()


in(x::AbstractVector, d::EmptyDomain) = false

# Arithemetic operations

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

in(x::AbstractVector, d::RnDomain) = true

# Arithemetic operations

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
end

Interval{T}(a::T = -1.0, b::T = 1.0) = Interval{T}(a, b)

Interval{T <: FloatingPoint}(::Type{T}) = Interval(zero(T), one(T))

in(x::AbstractVector, d::Interval) = in(x[1], d.a, d.b)

left(d::Interval) = d.a
right(d::Interval) = d.b

# Arithemetic operations

(+)(d::Interval, x::Number) = Interval(d.a+x, d.b+x)
(+)(x::Number, d::Interval) = d+x

(*)(d::Interval, x::Number) = Interval(x*d.a, x*d.b)
(*)(x::Number, d::Interval) = d*x

(/)(d::Interval, x::Number) = d * (1/x)

(==)(d1::Interval,d2::Interval) = (d1.a == d2.a) && (d1.b == d2.b)


show(io::IO, d::Interval) = print(io, "the interval [", d.a, ", ", d.b, "]")

const unitinterval = Interval()


###############################################################################################
### A circle
###############################################################################################

immutable Circle{T} <: AbstractDomain{2,T}
  radius    ::  T
  center    ::  Vector{T}

  Circle(radius::T = one(T), center::Vector{T} = zeros(T,2)) = new(radius, center) 
end

Circle{T}(radius::T = 1.0, center::Vector{T} = zeros(T,2)) = Circle{T}(radius,center)

in{T}(x::AbstractVector, c::Circle{T}) = (x[1]-c.center[1])^2 + (x[2]-c.center[2])^2 <= c.radius^2+10eps(T)

## Arithmetic operations

(+)(c::Circle, x::AbstractVector) = Circle(c.radius, c.center+x)
(+)(x::AbstractVector, c::Circle) = c+x

(*)(c::Circle, x::Number) = Circle(c.radius*x, c.center*x)
(*)(x::Number, c::Circle) = c*x

(/)(c::Circle, x::Number) = c * (1/x)

(==)(c1::Circle,c2::Circle) = (c1.radius == c2.radius) && (c1.center == c2.center)

show(io::IO, c::Circle) = print(io, "a circle of radius ", c.radius, " centered at ", c.center)

const unitcircle = Circle()



###############################################################################################
### A sphere
###############################################################################################

immutable Sphere{T} <: AbstractDomain{3,T}
  radius    ::  T
  center    ::  Vector{T}
end

Sphere{T}(radius::T = 1.0, center::Vector{T} = zeros(T,3)) = Sphere{T}(radius, center)

in(x::AbstractVector, s::Sphere) = (x[1]-s.center[1])^2 + (x[2]-s.center[2])^2 + (x[3]-s.center[3])^2 <= s.radius^2

## Arithmetic operations

(+)(s::Sphere, x::AbstractVector) = Sphere(s.radius, s.center+x)
(+)(x::AbstractVector, s::Sphere) = Sphere(s.radius, s.center+x)

(*)(s::Sphere, x::Number) = Sphere(s.radius * x, s.center * x)
(*)(x::Number, s::Sphere) = s*x

(/)(s::Sphere, x::Number) = s * (1/x)

(==)(s1::Sphere, s2::Sphere) = (s1.radius == s2.radius) && (s1.center == s2.center)

show(io::IO, s::Sphere) = print(io, "a sphere of radius ", s.radius, " centered at ", s.center)

const unitsphere = Sphere()


###############################################################################################
### An n-dimensional cube
###############################################################################################

immutable Cube{N,T} <: AbstractDomain{N,T}
  verts ::  Array{T,2}
end


Cube{T <: Number}(a::T, b::T) = Interval(a,b)

Cube{N,T}(left::NTuple{N,T}, right::NTuple{N,T}) = Cube{N,T}([[left...] [right...]])

Cube{T <: FloatingPoint}(::Type{T}) = Cube( (-one(T),-one(T),-one(T)), (one(T), one(T), one(T)))

Cube{T <: FloatingPoint}(n::Int, ::Type{T}) = Cube( tuple([-one(T) for i=1:n]...), tuple([one(T) for i=1:n]...))

Cube() = Cube(Float64)

Cube(n::Int) = Cube(n, Float64)

in{N}(x::AbstractVector, c::Cube{N}) = reduce(&, [in(x[j], c.verts[j,1], c.verts[j,2]) for j=1:N])

## Arithmetic operations

(+)(c::Cube, x::AbstractVector) = Cube(c.verts .+ x)
(+)(x::AbstractVector, c::Cube) = c+x

(*)(c::Cube, x::Number) = Cube(c.verts * x)
(*)(x::Number, c::Cube) = c*x

(/)(c::Cube, x::Number) = c * (1/x)

(==)(c1::Cube, c2::Cube) = (c1.verts == c2.verts)

show(io::IO, c::Cube{2}) = print(io, "the rectangle [", c.verts[1,1], ",", c.verts[1,2], "] x [", c.verts[2,1], ",", c.verts[2,2], "]")

show(io::IO, c::Cube{3}) = print(io, "the cube [", c.verts[1,1], ",", c.verts[1,2], "] x [", c.verts[2,1], ",", c.verts[2,2], "] x [", c.verts[3,1], ",", c.verts[3,2], "]")

show{N}(io::IO, c::Cube{N}) = print(io, "a ", N, "-dimensional cube")

const unitsquare = Cube(2)
const unitcube = Cube(3)



###############################################################################################
### A cylinder
###############################################################################################


immutable Cylinder{T} <: AbstractDomain{3,T}
  radius    ::  T
  length    ::  T
end

Cylinder{T}(radius::T = one(T), length::T = one(T)) = Cylinder{T}(radius, length)

in{T}(x::AbstractVector, c::Cylinder{T}) = in(x[1], zero(T), c.length) && (x[2]^2+x[3]^2 <= c.radius^2)

(==)(c1::Cylinder, c2::Cylinder) = (c1.radius == c2.radius) && (c1.length == c2.length)

show(io::IO, c::Cylinder) = print(io, "a cylinder of radius ", c.radius, " and length ", c.length)


###############################################################################################
### The union of two domains
###############################################################################################

# Type parameters N T D1 D2: dimension, numeric type, type of domain 1, type of domain 2.
# A stricter definition would use triangular dispatch:
# immutable DomainUnion{N, T, D1 <: AbstractDomain{N,T}, D2 <: AbstractDomain{N,T}} <: AbstractDomain{N,T}
immutable DomainUnion{N,T,D1,D2} <: AbstractDomain{N,T}
  d1    ::  D1
  d2    ::  D2
end

DomainUnion{N,T}(d1::AbstractDomain{N,T}, d2::AbstractDomain{N,T}) = DomainUnion{N,T,typeof(d1),typeof(d2)}(d1, d2)

join(d1::AbstractDomain, d2::AbstractDomain) = (d1 == d2 ? d1 : DomainUnion(d1,d2))

function join(d1::Interval, d2::Interval)
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
in(x::AbstractVector, d::DomainUnion) = in(x, d.d1) || in(x, d.d2)

function in(g::AbstractGrid, d::DomainUnion)
  z1 = in(g, d.d1)
  z2 = in(g, d.d2)
  z1 | z2
end

(+)(d1::AbstractDomain, d2::AbstractDomain) = join(d1,d2)
(|)(d1::AbstractDomain, d2::AbstractDomain) = join(d1,d2)

(==)(d1::DomainUnion, d2::DomainUnion) = (d1.d1 == d2.d1) && (d1.d2 == d2.d2)


function show(io::IO, d::DomainUnion)
    print(io, "A union of two domains: \n")
    print(io, "First domain: ", d.d1, "\n")
    print(io, "Second domain: ", d.d1, "\n")
end


###############################################################################################
### The intersection of two domains
###############################################################################################

immutable DomainIntersection{N,T,D1,D2} <: AbstractDomain{N,T}
  d1    ::  D1
  d2    ::  D2
end

DomainIntersection{N,T}(d1::AbstractDomain{N,T},d2::AbstractDomain{N,T}) = DomainIntersection{N,T,typeof(d1),typeof(d2)}(d1, d2)

# The intersection of two domains corresponds to a logical AND of their characteristic functions
in(x::AbstractVector, d::DomainIntersection) = in(x, d.d1) && in(x, d.d2)

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

(==)(d1::DomainIntersection, d2::DomainIntersection) = (d1.d1 == d2.d1) && (d1.d2 == d2.d2)

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
end

DomainDifference{N,T}(d1::AbstractDomain{N,T}, d2::AbstractDomain{N,T}) = DomainDifference{N,T,typeof(d1),typeof(d2)}(d1,d2)

# The difference between two domains corresponds to a logical AND NOT of their characteristic functions
in(x::AbstractVector, d::DomainDifference) = in(x, d.d1) && (~in(x, d.d2))

function in(g::AbstractGrid, d::DomainDifference)
    z1 = in(g, d.d1)
    z2 = in(g, d.d2)
    z1 & (~z2)
end

(-)(d1::AbstractDomain, d2::AbstractDomain) = DomainDifference(d1,d2)
(\)(d1::AbstractDomain, d2::AbstractDomain) = DomainDifference(d1,d2)

(==)(d1::DomainDifference, d2::DomainDifference) = (d1.d1 == d2.d1) && (d1.d2 == d2.d2)

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

function in(x::AbstractVector, d::RevolvedDomain)
    r = sqrt(x[2]^2+x[3])
    phi = atan2(x[2]/x[1])
    theta = acos(x[3]/r)
    in((x[1],r), d.d)
end

(==)(d1::RevolvedDomain, d2::RevolvedDomain) = (d1.d == d2.d)


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

in(x::AbstractVector, d::RotatedDomain) = in(d.rotationmatrix*x, d.d)

(==)(d1::RotatedDomain, d2::RotatedDomain) = (d1.d == d2.d) && (d1.angle == d2.angle) #&& (d1.rotationmatrix == d2.rotationmatrix)


###############################################################################################
### A scaled domain
###############################################################################################

# Note that the Euclidean plane is scaled, not just the domain itself.
# So the location of the origin matters. Two times a circle of radius 1 at a distance d of the origin
# becomes a circle of radius 2 at a distance 2d of the origin.
immutable ScaledDomain{N,T,D} <: AbstractDomain{N,T}
    d           ::  D
    scalefactor ::  T
end

ScaledDomain{N,T}(d::AbstractDomain{N,T}, scalefactor) = ScaledDomain{N,T,typeof(d)}(d, scalefactor)

function in(x::AbstractVector, d::ScaledDomain)
  in(x/d.scalefactor, d.d)
end

(*){N,T <: Number}(a::T, d::AbstractDomain{N,T}) = ScaledDomain(d,a)
(*){N,T <: Number}(d::AbstractDomain{N,T}, a::T) = a*d


###############################################################################################
### A translated domain
###############################################################################################

immutable TranslatedDomain{N,T,D} <: AbstractDomain{N,T}
    d       ::  D
    trans   ::  Vector{T}
end

TranslatedDomain{N,T}(d::AbstractDomain{N,T}, trans::Vector{T}) = TranslatedDomain{N,T,typeof(d)}(d, trans)

function in(x::AbstractVector, d::TranslatedDomain)
    in(x-d.trans, d.d)
end

(+){N,T}(d::AbstractDomain{N,T}, trans::AbstractVector{T}) = TranslatedDomain(d,trans)
(+){N,T}(trans::AbstractVector{T}, d::AbstractDomain{N,T}) = d + a


###############################################################################################
### A collection of domains
###############################################################################################

type DomainCollection{N,T} <: AbstractDomain{N,T}
    list    ::  Vector{AbstractDomain{N,T}}
end

DomainCollection{N,T}(d::AbstractDomain{N,T}) = DomainCollection{N,T}([d])

function in(x::AbstractVector, d::DomainCollection)
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

show(io::IO, d::DomainCollection) = print(io, "a collection of ", length(d.list), " domains")


function randomcircles(n)
    list = [Circle(0.2, (2*rand(2)-1)*0.8) for i=1:n]
    DC = DomainCollection(list[1])
    for i = 2:n
        push!(DC.list, list[i])
    end
    DC.box = FBox(-1.0, 1.0, -1.0, 1.0)
    DC
end


###
# The atomium
###

function atomium()
  sphere1 = Sphere(0.25)
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









