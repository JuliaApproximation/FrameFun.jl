abstract type AbstractSpace{N,T} end

basis(space::AbstractSpace) = space.basis
ndims{N,T}(::Type{AbstractSpace{N,T}}) = N
ndims{S <: AbstractSpace}(::Type{S}) = ndims(supertype(S))
ndims{N,T}(::AbstractSpace{N,T}) = N
numtype(s::AbstractSpace) = real(eltype(s))
eltype{N,T}(::Type{AbstractSpace{N,T}}) = T
eltype{B <: AbstractSpace}(::Type{B}) = eltype(supertype(B))
FunctionSet{N,T}(space::AbstractSpace{N,T},n) = resize(basis(space),n)

struct FunctionSpace{N,T} <: AbstractSpace{N,T}
  basis   ::    FunctionSet{N,T}

  FunctionSpace{N,T}(basis::FunctionSet) where {N,T} = new(basis)
  FunctionSpace{N,T}(basis::FunctionSet, dom::Domain) where {N,T} = new(rescale(basis,leftendpoint(dom),rightendpoint(dom)))
end

FunctionSpace{N,T}(basis::FunctionSet{N,T}) = FunctionSpace{N,T}(basis)
FunctionSpace{N,T}(basis::FunctionSet{N,T}, dom::Domain) = FunctionSpace{N,T}(basis,dom)
# place somewhere else?
FourierSpace(left=0,right=1) = FunctionSpace(FourierBasis(0),interval(left,right))
ChebyshevSpace(left=-1,right=1) = FunctionSpace(ChebyshevBasis(0),interval(left,right))
# place somewhere else?
boundingbox{N,T}(f::FunctionSet{N,T}) = boundingbox(left(f),right(f))

boundingbox(space::FunctionSpace) = boundingbox(space.basis)

name(space::AbstractSpace) = "Space of "*name(FunctionSet(space,0))
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

Base.promote_eltype{N,T}(space::AbstractSpace{N,T}, ::Type{T}) = space
Base.promote_eltype{N,T1, T2}(space::AbstractSpace{N,T1}, ::Type{T2}) = FunctionSpace(promote_eltype(basis(space), T2))
Base.promote{N,T}(space1::AbstractSpace{N,T},space2::AbstractSpace{N,T}) = (space1,space2)
function Base.promote{N,T1,T2}(space1::FunctionSpace{N,T1}, space2::FunctionSpace{N,T2})
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
