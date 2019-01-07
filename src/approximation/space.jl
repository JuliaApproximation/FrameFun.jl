abstract type AbstractSpace{S,T}
end

basis(space::AbstractSpace) = space.basis

dimension(space::B) where{B<:AbstractSpace} = dimension(basis(space))
domaintype(space::B) where{B<:AbstractSpace} = domaintype(basis(space))
codomaintype(space::B) where{B<:AbstractSpace} = codomaintype(basis(space))

Dictionary(space::AbstractSpace{S,T},n) where {S,T}= resize(basis(space),n) 

struct GenericFunctionSpace{S,T} <: AbstractSpace{S,T}
    basis   ::    Dictionary{S,T}

    GenericFunctionSpace{S,T}(basis::Dictionary) where {S,T} = new(basis)
end

GenericFunctionSpace(basis::Dictionary{S,T}) where {S,T} = GenericFunctionSpace{S,T}(basis)
GenericFunctionSpace(basis::Dictionary, dom::Domain) =
    GenericFunctionSpace(rescale(basis, boundingbox(dom)))

# place somewhere else?
FourierSpace(left::Real=0,right::Real=1) =
    GenericFunctionSpace(FourierBasis(0), Interval(left,right))
ChebyshevSpace(left::Real=-1,right::Real=1) =
    GenericFunctionSpace(ChebyshevBasis(0), Interval(left,right))
# place somewhere else?
boundingbox(f::Dictionary) = boundingbox(support(f))

boundingbox(space::GenericFunctionSpace) = boundingbox(space.basis)

name(space::AbstractSpace) = "Space of "*name(Dictionary(space,0))

"Tensorproduct of function space"
tensorproduct(space::GenericFunctionSpace) = space
tensorproduct(space::GenericFunctionSpace, n::Int) = tensorproduct([space for i=1:n]...)
tensorproduct(space1::GenericFunctionSpace, space2::GenericFunctionSpace, spaces::GenericFunctionSpace...) =
#tensorproduct of functionspace is tensorproduct of bases and cartesian product of boundingboxes
    tensorproduct(GenericFunctionSpace(tensorproduct(basis(space1),basis(space2)), cartesianproduct(boundingbox(space1),boundingbox(space2))), spaces...)

"""
Addition of multiple function spaces

The intervals of the function spaces are transformed to the union of
the bounding boxes of the function spaces.
"""
add(space::GenericFunctionSpace) = space
add(space::GenericFunctionSpace, n::Int) = ⊕([space for i=1:n]...)
add(space1::GenericFunctionSpace, space2::GenericFunctionSpace, spaces::GenericFunctionSpace...) =
    add(GenericFunctionSpace(basis(space1)⊕basis(space2),union(boundingbox(space1),boundingbox(space2))), spaces...)
⊕(args::GenericFunctionSpace...) = add(args...)

promote_domaintype(space::AbstractSpace{T1}, ::Type{T2}) where {T1,T2} = GenericFunctionSpace(promote_domaintype(basis(space), T2))

# for op in (:left, :right)
#   @eval begin
#     $op(f::GenericFunctionSpace) = $op(boundingbox(f))
#   end
# end

# for op in (:is_basis, :is_frame, :isorthogonal, :is_biorthogonal,
#   :has_derivative, :has_antiderivative, :has_interpolationgrid, :has_transform, :has_extension)
#   @eval begin
#     $op(f::GenericFunctionSpace) = $op(basis(f))
#   end
# end
