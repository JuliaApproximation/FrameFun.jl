# box.jl

"An BBox is an N-dimensional box specified by its bottom-left and top-right vertices."
immutable BBox{N,T}
    left        ::  Vec{N,T}
    right       ::  Vec{N,T}
end

dim{N,T}(::Type{BBox{N,T}}) = N
dim(b::BBox) = dim(typeof(b))

eltype{N,T}(::Type{BBox{N,T}}) = T


typealias BBox1{T} BBox{1,T}
typealias BBox2{T} BBox{2,T}
typealias BBox3{T} BBox{3,T}
typealias BBox4{T} BBox{4,T}

# Dimension-specific constructors
BBox{T}(a::T, b::T) = BBox(Vec{1,T}(a), Vec{1,T}(b))
BBox{T}(a::T, b::T, c::T, d::T) = BBox(Vec{2,T}(a,c), Vec{2,T}(b,d))
BBox{T}(a::T, b::T, c::T, d::T, e::T, f::T) = BBox(Vec{3,T}(a,c,e), Vec{3,T}(b,d,f))

BBox{N,T}(left::NTuple{N,T}, right::NTuple{N,T}) = BBox{N,T}(Vec{N,T}(left), Vec{N,T}(right))

BBox(left::AbstractVector, right::AbstractVector) = BBox(left, right, Val{length(left)})
BBox{N,T}(left::AbstractVector{T}, right::AbstractVector{T}, ::Type{Val{N}}) =
    BBox{N,T}(Vec{N,T}(left), Vec{N,T}(right))

left(b::BBox) = b.left
left(b::BBox1) = b.left[1]
left(b::BBox, dim) = b.left[dim]

right(b::BBox) = b.right
right(b::BBox1) = b.right[1]
right(b::BBox, dim) = b.right[dim]

size(b::BBox, dim) = right(b, dim) - left(b, dim)


# Extend a box by a factor of t[i] in each dimension
function extend{N,T}(b::BBox{N,T}, t::Vec{N,T})
    r = Vec{N,T}( [ t[i]*size(b,i) for i in 1:N ] )
    BBox(left(b), left(b) + r)
end

in(x, b::BBox) = reduce(&, [in(x,b,i) for i=1:dim(b)])
in{N,T}(x, b::BBox{N,T}, dim) = (x >= left(b,dim)-10eps(T)) && (x <= right(b,dim)+10eps(T))

## Arithmetic operations

(*)(a::Number,b::BBox) = BBox(a*b.left, a*b.right)
(*)(b::BBox, a::Number) = a*b
(/)(b::BBox, a::Number) = BBox(b.left/a, b.right/a)
(+)(b::BBox, a::Vector) = BBox(b.left+a, b.right+a)
(-)(b::BBox, a::Vector) = BBox(b.left-a, b.right-a)

# Logical operations on boxes: union and intersection
(|)(b1::BBox, b2::BBox) = BBox(min(left(b1),left(b2)), max(right(b1),right(b2)))
(+)(b1::BBox, b2::BBox) = b1 | b2
(&)(b1::BBox, b2::BBox) = BBox(max(left(b1),left(b2)), min(right(b1),right(b2)))

join(b1::BBox, b2::BBox) = b1 | b2
intersect(b1::BBox, b2::BBox) = b1 & b2

(==)(b1::BBox, b2::BBox) = (b1.left == b2.left) && (b1.right == b2.right)

(≈)(b1::BBox, b2::BBox) = (b1.left ≈ b2.left) && (b1.right ≈ b2.right)


show(io::IO, c::BBox{1}) = print(io, "the interval [", left(c, 1), ",", right(c, 1), "]")

show(io::IO, c::BBox{2}) = print(io, "the rectangular box [", left(c, 1), ",", right(c, 1), "] x [", left(c, 2), ",", right(c, 2), "]")

show(io::IO, c::BBox{3}) = print(io, "the cube [", left(c, 1), ",", right(c, 1), "] x [", left(c, 2), ",", right(c, 2), "] x [", left(c, 3), ",", right(c, 23), "]")

# Define the unit box
const unitbox1 = BBox(-1.0, 1.0)
const unitbox2 = BBox(-1.0, 1.0, -1.0, 1.0)
const unitbox3 = BBox(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)


