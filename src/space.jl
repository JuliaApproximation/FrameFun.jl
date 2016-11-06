abstract AbstractSpace{N,T}

basis(space::AbstractSpace) = space.basis
ndims{N,T}(::Type{AbstractSpace{N,T}}) = N
ndims{N,T}(::AbstractSpace{N,T}) = N
numtype(s::AbstractSpace) = real(eltype(s))
eltype{N,T}(::Type{AbstractSpace{N,T}}) = T
eltype{B <: AbstractSpace}(::Type{B}) = eltype(supertype(B))
FunctionSet{N,T}(space::AbstractSpace{N,T},n) = resize(basis(space),n)

immutable FunctionSpace{N,T,RT} <: AbstractSpace{N,T}
  basis   ::    FunctionSet{N,T}
  bbox    ::    BBox{N,RT}

  FunctionSpace(basis::FunctionSet) = new(basis,boundingbox(basis))
  FunctionSpace(basis::FunctionSet, bbox::BBox) = new(rescale(basis,left(bbox),right(bbox)),bbox)
end

FunctionSpace{N,T}(basis::FunctionSet{N,T}) = FunctionSpace{N,T,real(T)}(basis)
FunctionSpace{N,T}(basis::FunctionSet{N,T}, bbox::BBox) = FunctionSpace{N,T,real(T)}(basis,bbox)
# place somewhere else?
FourierSpace(left=0,right=1) = FunctionSpace(FourierBasis(0),BBox(left,right))
ChebyshevSpace(left=-1,right=1) = FunctionSpace(ChebyshevBasis(0),BBox(left,right))
# place somewhere else?
boundingbox{N,T}(f::FunctionSet{N,T}) = BBox{N,real(T)}(left(f),right(f))

boundingbox(space::FunctionSpace) = space.bbox

"Tensorproduct of function space"
tensorproduct(space::FunctionSpace) = space
tensorproduct(space::FunctionSpace, n::Int) = tensorproduct([space for i=1:n]...)
tensorproduct(space1::FunctionSpace, space2::FunctionSpace, spaces::FunctionSpace...) =
#tensorproduct of functionspace is tensorproduct of bases and boundingboxes
  tensorproduct(FunctionSpace(tensorproduct(basis(space1),basis(space2)), tensorproduct(boundingbox(space1),boundingbox(space2))), spaces...)

"""
Addition of multiple function spaces

The intervals of the function spaces are transformed to the union of
the bounding boxes of the function spaces.
"""
add(space::FunctionSpace) = space
add(space::FunctionSpace, n::Int) = ⊕([space for i=1:n]...)
add(space1::FunctionSpace, space2::FunctionSpace, spaces::FunctionSpace...) =
  add(FunctionSpace(basis(space1)⊕basis(space2),union(boundingbox(space1),boundingbox(space2))), spaces...)
⊕(args::FunctionSpace...) = add(args...)

promote_eltype{N,T}(space::AbstractSpace{N,T}, ::Type{T}) = s
promote{N,T}(space1::AbstractSpace{N,T},space2::AbstractSpace{N,T}) = (space1,space2)
function promote{N,T1,T2}(set1::FunctionSpace{N,T1}, set2::FunctionSpace{N,T2})
  T = promote_type(T1,T2)
  (promote_eltype(space1,T), promote_eltype(space2,T))
end

# for op in (:left, :right)
#   @eval begin
#     $op(f::FunctionSpace) = $op(boundingbox(f))
#   end
# end

# for op in (:is_basis, :is_frame, :is_orthogonal, :is_biorthogonal,
#   :has_derivative, :has_antiderivative, :has_grid, :has_transform, :has_extension)
#   @eval begin
#     $op(f::FunctionSpace) = $op(basis(f))
#   end
# end

function test_contructor()
  println("test 1")
  f = FourierBasis(121)⊗FourierBasis(121)
  fs = FunctionSpace(f)
  @assert left(fs)==Vector([0,0])
  @assert right(fs)==Vector([1,1])

  println("test 2")
  f = FourierBasis(121)
  fs = FunctionSpace(f)
  @assert left(fs)==0
  @assert right(fs)==1

  println("test 3")
  f = FourierBasis(121)⊗FourierBasis(121)
  b = BBox(-2.,0,-1,1)
  fs = FunctionSpace(f,b)
  @assert left(fs)==Vector([-2.,-1])
  @assert right(fs)==Vector([0,1])
end
