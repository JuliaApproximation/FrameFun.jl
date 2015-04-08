# box.jl

# The definition of a box is a 2D-array with the bottom-left and top-right vertices as columns.
immutable FBox{N, T <: FloatingPoint}
	vertices::Array{T,2}
end

dim{N,T}(::FBox{N,T}) = N
dim{N,T}(::Type{FBox{N,T}}) = N
dim{B <: FBox}(::Type{B}) = dim(super(B))

numtype{N,T}(::FBox{N,T}) = T
numtype{N,T}(::Type{FBox{N,T}}) = T
numtype{B <: FBox}(::Type{B}) = numtype(super(B))

typealias FBox1{T <: FloatingPoint} FBox{1,T}
typealias FBox2{T <: FloatingPoint} FBox{2,T}
typealias FBox3{T <: FloatingPoint} FBox{3,T}
typealias FBox4{T <: FloatingPoint} FBox{4,T}

# A type-safe general constructor is hard: how can N and T be inferred?
FBox{T}(a::Vector{T}, b::Vector{T}) = FBox{length(a),T}([a b])

# Dimension-specific constructors
FBox{T}(a::T, b::T) = FBox{1,T}([a b])
FBox{T}(a::T, b::T, c::T, d::T) = FBox{2,T}([a b; c d])
FBox{T}(intervalx::Vector{T}) = FBox{1,T}([intervalx[1] intervalx[2]])
FBox{T}(intervalx::Vector{T}, intervaly::Vector{T}) = FBox{2,T}([intervalx[1] intervalx[2]; intervaly[1] intervaly[2]])
FBox{T}(intervalx::Vector{T}, intervaly::Vector{T}, intervalz::Vector{T}) = FBox{3,T}([intervalx[1] intervalx[2]; intervaly[1] intervaly[2]; intervalz[1] intervalz[2]])

emptybox(N,T) = FBox{N,T}(zeros(T,N,2))

# Other constructors are specific for each dimension
#FBox{N,T}(a::Vector{T},b::Vector{T}) = FBox{N,T}([a b])
#FBox{T,1}(a::T, b::T) = FBox{T,1}([a b])


left(b::FBox) = b.vertices[:,1]
left(b::FBox, dim) = b.vertices[dim,1]
right(b::FBox) = b.vertices[:,2]
right(b::FBox, dim) = b.vertices[dim,2]

size(b::FBox, dim) = b.vertices[dim,2]-b.vertices[dim,1]

vertices(b::FBox) = b.vertices

# Extend a box by a factor of t[i] in each dimension
function extend{N,T}(b::FBox{N,T}, t::NTuple{N,T})
	vertices = copy(b.vertices)
	for i=1:size(vertices,1)
		vertices[i,2] = vertices[i,1] + t[i]*(vertices[i,2]-vertices[i,1])
	end
	FBox{N,T}(vertices)
end

in(x, b::FBox) = reduce(&, [in(x,b,i) for i=1:dim(b)])
in{N,T}(x, b::FBox{N,T}, dim) = (x >= left(b,dim)-10eps(T)) && (x <= right(b,dim)+10eps(T))

## Arithemetic operations

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


