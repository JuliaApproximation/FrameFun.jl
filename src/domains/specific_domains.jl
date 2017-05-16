# specific_domains.jl

################################################################################
### An empty domain
################################################################################

immutable EmptyDomain{N} <: AbstractDomain{N}
end

EmptyDomain() = EmptyDomain{1}()
EmptyDomain(n::Int) = EmptyDomain{n}()
# The typesafe constructor is EmptyDomain{N}()

indomain(x, d::EmptyDomain) = false

# Arithmetic operations

(+)(d::EmptyDomain, x::Number) = d

(*)(d::EmptyDomain, x::Number) = d


show(io::IO, d::EmptyDomain) = print(io, "an empty domain")


################################################################################
### The space R^n
################################################################################

immutable RnDomain{N} <: AbstractDomain{N}
end

RnDomain() = RnDomain{1}()
RnDomain(n::Int) = RnDomain{n}()
# The typesafe constructor is RnDomain{N}()

indomain(x, d::RnDomain) = true

# Arithmetic operations

(+)(d::RnDomain, x::Number) = d

(*)(d::RnDomain, x::Number) = d


show(io::IO, e::RnDomain) = print(io, "the ", ndims(e), "-dimensional Euclidean space")



################################################################################
### An interval
################################################################################

immutable Interval{T} <: AbstractDomain{1}
    a     ::  T
    b     ::  T

    function Interval(a = -1, b = 1)
      @assert a<b
      new(a,b)
    end
end

Interval() = Interval{Int}()

Interval{T}(::Type{T}) = Interval{T}()

Interval{T <: Number}(a::T, b::T) = Interval{T}(a, b)
Interval{S <: Number, T <: Number}(a::S, b::T) = Interval(promote(a,b)...)


indomain(x, d::Interval) = in(x, d.a, d.b)

left(d::Interval) = d.a
right(d::Interval) = d.b

length(d::Interval) = right(d)-left(d)

# Arithmetic operations

(+)(d::Interval, x::Number) = Interval(d.a+x, d.b+x)

(*)(a::Number, d::Interval) = Interval(a*d.a, a*d.b)
(*)(d::Interval, a::Number) = a * d


boundingbox(d::Interval) = BBox(left(d), right(d))

show(io::IO, d::Interval) = print(io, "the interval [", d.a, ", ", d.b, "]")

const unitinterval = Interval()

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

function intersect(d1::Interval, d2::Interval)
    a = left(d1)
    b = right(d1)
    c = left(d2)
    d = right(d2)

    if (b < c) || (a > d)
        EmptyDomain(Val{1}())
    else
        Interval(max(a, c), min(b, d))
    end
end


################################################################################
### The unit ball in N dimensions
################################################################################

immutable UnitBall{N} <: AbstractDomain{N}
end

indomain(x, ::UnitBall{1}) = -1 <= x <= 1
indomain(x, ::UnitBall{2}) = x[1]^2+x[2]^2 <= 1
indomain(x, ::UnitBall{3}) = x[1]^2+x[2]^2+x[3]^2 <= 1

indomain{N}(x, ::UnitBall{N}) = sum(map(t->t^2, x)) <= 1

boundingbox{N}(::UnitBall{N}) = BBox{N,Int}(-ones(SVector{N,Int}), ones(SVector{N,Int}))

Disk() = UnitBall{2}()
Disk(radius) = radius * Disk()
Disk(radius, center) = radius * Disk() + center

show(io::IO, d::UnitBall) = print(io, "the $(ndims(d))-dimensional unit ball")


# ################################################################################
# ### A disk
# ################################################################################
#
# immutable Disk{S,T} <: AbstractDomain{2}
#     radius    ::  S
#     center    ::  SVector{2,T}
#
#     Disk(radius = one(S), center = SVector(0,0)) = new(radius, center)
# end
#
# Disk() = Disk{Int,Float64}()
# Disk{T}(::Type{T}) = Disk{T,T}()
#
# Disk{T}(radius::T) = Disk{T,T}(radius)
# Disk{S,T}(radius::S, center::SVector{2,T}) = Disk{S,T}(radius, center)
# Disk(radius, center::AbstractVector) = Disk(radius, SVector{2}(center))
#
# indomain(x, c::Disk) = (x[1]-c.center[1])^2 + (x[2]-c.center[2])^2 <= c.radius^2
#
# ## Arithmetic operations
#
# (+)(c::Disk, x::SVector{2}) = Disk(c.radius, c.center+x)
#
# (*)(c::Disk, x::Number) = Disk(c.radius*x, c.center*x)
#
#
# boundingbox(c::Disk) = BBox((c.center[1]-c.radius,c.center[2]-c.radius),(c.center[1]+c.radius,c.center[2]+c.radius))
#
# show(io::IO, c::Disk) = print(io, "a disk of radius ", c.radius, " centered at ", c.center)

const unitdisk = Disk()



################################################################################
### A 3D ball
################################################################################

immutable Ball{S,T} <: AbstractDomain{3}
    radius    ::  S
    center    ::  SVector{3,T}

    Ball(radius = one(S), center = zeros(SVector{3,T})) = new(radius, center)
end

Ball() = Ball{Int,Float64}()
Ball{T}(::Type{T}) = Ball{T,T}()

Ball{T}(radius::T) = Ball{T,T}(radius)
Ball{S,T}(radius::S, center::SVector{3,T}) = Ball{S,T}(radius, center)
Ball(radius, center::AbstractVector) = Ball(radius, SVector{3}(center))


indomain(x, s::Ball) = (x[1]-s.center[1])^2 + (x[2]-s.center[2])^2 + (x[3]-s.center[3])^2 <= s.radius^2

## Arithmetic operations

(+)(s::Ball, x::SVector{3}) = Ball(s.radius, s.center+x)

(*)(s::Ball, x::Number) = Ball(s.radius * x, s.center * x)


boundingbox(c::Ball) = BBox((c.center[1]-c.radius,c.center[2]-c.radius,c.center[3]-c.radius),(c.center[1]+c.radius,c.center[2]+c.radius,c.center[3]+c.radius))

show(io::IO, s::Ball) = print(io, "a ball of radius ", s.radius, " centered at ", s.center)

const unitball = Ball()





################################################################################
### An n-dimensional cube
################################################################################

Cube() = Cube(Val{1})

Cube(::Type{Val{1}}) = Interval()
Cube{N}(::Type{Val{N}}) = Interval() ⊗ Cube(Val{N-1})
Cube(n::Int) = Cube(Val{n})

Cube(left::Number, right::Number) = Interval(left, right)
Cube(left, right) = tensorproduct(map(Interval, left, right)...)

rectangle(a, b, c, d) = Interval(a,b) ⊗ Interval(c,d)

cube(a, b, c, d, e, f) = Interval(a,b) ⊗ Interval(c,d) ⊗ Interval(e,f)


const unitsquare = Cube(Val{2})
const unitcube = Cube(Val{3})



################################################################################
### A cylinder
################################################################################


cylinder(radius = 1, length = 1) = Disk(radius) ⊗ Interval(0,length)



##########################################################################
### Assorted Domains
##########################################################################


function randomcircles(n, radius = 0.3)
    list = AbstractDomain[Disk(radius, (2*rand(2)-1)*0.8) for i=1:n]
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
