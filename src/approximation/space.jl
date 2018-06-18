abstract type AbstractSpace{S,T}
end

basis(space::AbstractSpace) = space.basis

dimension(space::B) where{B<:AbstractSpace} = dimension(basis(space))
domaintype(space::B) where{B<:AbstractSpace} = domaintype(basis(space))
codomaintype(space::B) where{B<:AbstractSpace} = codomaintype(basis(space))

Dictionary{S,T}(space::AbstractSpace{S,T},n) = resize(basis(space),n)

struct FunctionSpace{S,T} <: AbstractSpace{S,T}
    basis   ::    Dictionary{S,T}

    FunctionSpace{S,T}(basis::Dictionary) where {S,T} = new(basis)
end

FunctionSpace(basis::Dictionary{S,T}) where {S,T} = FunctionSpace{S,T}(basis)
FunctionSpace(basis::Dictionary, dom::Domain) =
    FunctionSpace(rescale(basis, infimum(dom), supremum(dom)))

# place somewhere else?
FourierSpace(left::Real=0,right::Real=1) =
    FunctionSpace(FourierBasis(0), interval(left,right))
ChebyshevSpace(left::Real=-1,right::Real=1) =
    FunctionSpace(ChebyshevBasis(0), interval(left,right))
# place somewhere else?
boundingbox(f::Dictionary) = boundingbox(infimum(support(f)), supremum(support(f)))

boundingbox(space::FunctionSpace) = boundingbox(space.basis)

name(space::AbstractSpace) = "Space of "*name(Dictionary(space,0))

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
