abstract type AbstractSpace{T}
end

basis(space::AbstractSpace) = space.basis

dimension(space::B) where{B<:AbstractSpace} = dimension(basis(space))
domaintype(space::B) where{B<:AbstractSpace} = domaintype(basis(space))
rangetype(space::B) where{B<:AbstractSpace} = rangetype(basis(space))

FunctionSet{T}(space::AbstractSpace{T},n) = resize(basis(space),n)

struct FunctionSpace{T} <: AbstractSpace{T}
    basis   ::    FunctionSet{T}

    FunctionSpace{T}(basis::FunctionSet) where {T} = new(basis)
    FunctionSpace{T}(basis::FunctionSet, dom::Domain) where {T} = new(rescale(basis,leftendpoint(dom),rightendpoint(dom)))
end

FunctionSpace(basis::FunctionSet{T}) where {T} = FunctionSpace{T}(basis)
FunctionSpace(basis::FunctionSet{T}, dom::Domain) where {T} = FunctionSpace{T}(basis,dom)
# place somewhere else?
FourierSpace(left=0,right=1) = FunctionSpace(FourierBasis(0), interval(left,right))
ChebyshevSpace(left=-1,right=1) = FunctionSpace(ChebyshevBasis(0), interval(left,right))
# place somewhere else?
boundingbox(f::FunctionSet{T}) where {T} = boundingbox(left(f), right(f))

boundingbox(space::FunctionSpace) = boundingbox(space.basis)

name(space::AbstractSpace) = "Space of "*name(FunctionSet(space,0))

"Tensorproduct of function space"
tensorproduct(space::FunctionSpace) = space
tensorproduct(space::FunctionSpace, n::Int) = tensorproduct([space for i=1:n]...)
tensorproduct(space1::FunctionSpace, space2::FunctionSpace, spaces::FunctionSpace...) =
#tensorproduct of functionspace is tensorproduct of bases and cartesian product of boundingboxes
    tensorproduct(FunctionSpace(tensorproduct(basis(space1),basis(space2)), cartesianproduct(boundingbox(space1),boundingbox(space2))), spaces...)

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

promote_domaintype(space::AbstractSpace{T1}, ::Type{T2}) where {T1,T2} = FunctionSpace(promote_domaintype(basis(space), T2))

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
