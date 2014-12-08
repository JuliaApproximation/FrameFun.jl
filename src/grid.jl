# grid.jl

immutable Grid{N,T <: FloatingPoint}
	ranges::NTuple{N,FloatRange{T}}
end

dim{N}(g::Grid{N}) = N
numtype{N,T}(g::Grid{N,T}) = T

typealias Grid1{T <: FloatingPoint} Grid{1,T}
typealias Grid2{T <: FloatingPoint} Grid{2,T}
typealias Grid3{T <: FloatingPoint} Grid{3,T}
typealias Grid4{T <: FloatingPoint} Grid{4,T}

## Constructors

# Annotate type of first variable in order to avoid ambiguity later on
Grid(xrange::FloatRange) = Grid((xrange,))
Grid(xrange::FloatRange,yrange) = Grid((xrange,yrange))
Grid(xrange::FloatRange,yrange,zrange) = Grid((xrange,yrange,zrange))
Grid(xrange::FloatRange,yrange,zrange,trange) = Grid((xrange,yrange,zrange,trange))

# Create a grid from a box using m[i] points in each dimension
Grid{N,T}(b::FBox{N,T}, m::NTuple{N}) = Grid{N,T}(tuple([linrange(left(b,i),right(b,i), m[i]) for i=1:N]...))

periodicgrid{N,T}(b::FBox{N,T}, m) = Grid{N,T}(tuple([linrange(b.vertices[i,1],b.vertices[i,2]-(b.vertices[i,2]-b.vertices[i,1])/m[i], m[i]) for i=1:N]...))

# Create a subgrid of the given grid, with the first n[i] points in each dimension
subgrid{N,T}(g::Grid{N,T}, n) = Grid{N,T}(tuple([linrange(left(g,i),range(g,i)[n[i]],n[i]) for i=1:N]...))



size(g::Grid, j) = length(g.ranges[j])
size{N}(g::Grid{N}) = tuple([length(g.ranges[j]) for j=1:N]...)

left(g::Grid, j) = start(g.ranges[j])
left{N}(g::Grid{N}) = [left(g, j) for j=1:N]

right(g::Grid, j) = last(g.ranges[j])
right{N}(g::Grid{N}) = [right(g, j) for j=1:N]

range(g::Grid, j) = g.ranges[j]
range(g::Grid) = g.ranges

# Create a box that encapsulates the grid
box{N,T}(g::Grid{N,T}) = FBox{N,T}(left(g), right(g))

getindex{N}(g::Grid{N}, i...) = [g.ranges[k][i[k]] for i=1:N]
# getindex(g::Grid1, i::Int) = g.ranges[1][i]		# special case that does not create a vector


## Arithmetic operations
(==)(g1::Grid,g2::Grid) = g1.ranges == g2.ranges


