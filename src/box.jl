# box.jl

"An FBox is an N-dimensional box specified by its bottom-left and top-right vertices."
immutable FBox{N,T}
	vertices   ::  Mat{N,2,T}
end

dim{N,T}(::FBox{N,T}) = N
dim{N,T}(::Type{FBox{N,T}}) = N
dim{B <: FBox}(::Type{B}) = dim(super(B))

numtype{N,T}(::FBox{N,T}) = T
numtype{N,T}(::Type{FBox{N,T}}) = T
numtype{B <: FBox}(::Type{B}) = numtype(super(B))

typealias FBox1{T} FBox{1,T}
typealias FBox2{T} FBox{2,T}
typealias FBox3{T} FBox{3,T}
typealias FBox4{T} FBox{4,T}

FBox{N,T}(a::Vec{N,T}, b::Vec{N,T}) = FBox(Mat{N,2,T}((a,b)))

# Dimension-specific constructors
FBox{T}(a::T, b::T) = FBox(Vec{1,T}(a), Vec{1,T}(b))
FBox{T}(a::T, b::T, c::T, d::T) = FBox(Vec{2,T}(a,c), Vec{2,T}(b,d))
FBox{T}(a::T, b::T, c::T, d::T, e::T, f::T) = FBox(Vec{3,T}(a,c,e), Vec{3,T}(b,d,f))

left{N,T}(b::FBox{N,T}) = Vec{N,T}([b.vertices[i,1] for i in 1:N])
left(b::FBox, dim) = b.vertices[dim,1]

right{N,T}(b::FBox{N,T}) = Vec{N,T}([b.vertices[i,2] for i in 1:N])
right(b::FBox, dim) = b.vertices[dim,2]

size(b::FBox, dim) = right(b, dim) - left(b, dim)

vertices(b::FBox) = b.vertices

# Extend a box by a factor of t[i] in each dimension
function extend{N,T}(b::FBox{N,T}, t::Vec{N,T})
    r = Vec{N,T}( [ t[i]*size(b,i) for i in 1:N ] )
    FBox(left(b), left(b) + r)
end

in(x, b::FBox) = reduce(&, [in(x,b,i) for i=1:dim(b)])
in{N,T}(x, b::FBox{N,T}, dim) = (x >= left(b,dim)-10eps(T)) && (x <= right(b,dim)+10eps(T))

## Arithmetic operations

(*)(a::Number,b::FBox) = FBox(a*b.vertices)
(*)(b::FBox, a::Number) = a*b
(/)(b::FBox, a::Number) = FBox(b.vertices/a)
(+)(b::FBox, a::Vector) = FBox(left(b)+a, right(b)+a)
(-)(b::FBox, a::Vector) = FBox(b.vertices .- a)

# Logical operations on boxes: union and intersection
(|)(b1::FBox, b2::FBox) = FBox(min(left(b1),left(b2)), max(right(b1),right(b2)))
(+)(b1::FBox, b2::FBox) = b1|b2
(&)(b1::FBox, b2::FBox) = FBox(max(left(b1),left(b2)), min(right(b1),right(b2)))

(==)(b1::FBox, b2::FBox) = b1.vertices == b2.vertices

show(io::IO, b::FBox) = print(io, "The box ", b.vertices)

# Define the unit box
const unitbox1 = FBox(-1.0, 1.0)
const unitbox2 = FBox([-1.0, -1.0], [1.0, 1.0])
const unitbox3 = FBox([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0])
const unitbox4 = FBox([-1.0, -1.0, -1.0, -1.0], [1.0, 1.0, 1.0, 1.0])


