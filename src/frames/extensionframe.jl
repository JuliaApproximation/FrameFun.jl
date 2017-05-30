# extensionframe.jl

"""
An ExtensionFrame is the restriction of a basis to a subset of its domain. This
results in a frame that implicitly represents extensions of functions on the
smaller set to the larger set.
"""
struct ExtensionFrame{N,T} <: DerivedSet{N,T}
    domain      ::  Domain{N}
    basis       ::  FunctionSet{N,T}

    function ExtensionFrame{N,T}(domain::Domain, basis::FunctionSet) where {N,T}
        @assert is_basis(basis)
        new(domain, basis)
    end
end

ExtensionFrame{N,T}(domain::Domain{N}, basis::FunctionSet{N,T}) =
    ExtensionFrame{N,T}(domain, basis)

# superset is the function for DerivedSet's to obtain the underlying set
superset(f::ExtensionFrame) = f.basis

basis(f::ExtensionFrame) = f.basis
domain(f::ExtensionFrame) = f.domain

similar_set(f::ExtensionFrame, set::FunctionSet) = ExtensionFrame(domain(f), set)

is_frame(f::ExtensionFrame) = true

# The following properties do not hold for extension frames
# - there is no interpolation grid
has_grid(f::ExtensionFrame) = false
# - there is no unitary transform
has_transform(f::ExtensionFrame) = false
has_transform(f::ExtensionFrame, dgs) = false
# - there is no antiderivative (in general)
has_antiderivative(f::ExtensionFrame) = false

name(f::ExtensionFrame) = "An extension frame of " * name(f.basis)

in_support(f::ExtensionFrame, idx::Int, x) = x ∈ domain(f)

is_compatible(d1::ExtensionFrame, d2::ExtensionFrame) = is_compatible(basis(d1),basis(d2))

function (*)(d1::ExtensionFrame, d2::ExtensionFrame, args...)
    @assert is_compatible(d1,d2)
    (mset, mcoef) = (*)(basis(d1),basis(d2),args...)
    df = ExtensionFrame(domain(d1) ∩ domain(d2),mset)
    (df, mcoef)
end

# Should we check whether x lies in the domain?
call_set_expansion(e::SetExpansion, s::ExtensionFrame, coef, x...) = eval_expansion(basis(s), coef, x...)

eval_element(s::ExtensionFrame, idx::Int, x) = eval_element(basis(s), idx, x)

grid(f::ExtensionFrame) = subgrid(grid(basis(f)),domain(f))


"""
Make an ExtensionFrame, but match tensor product domains with tensor product sets
in a suitable way.

For example: an interval ⊗ a disk (= a cylinder) combined with a 3D Fourier series, leads to a
tensor product of a Fourier series on the interval ⊗ a 2D Fourier series on the disk.
"""
function extensionframe(domain::ProductDomain, basis::TensorProductSet)
    ExtensionFrames = FunctionSet[]
    dc = 1
    for i = 1:nb_elements(domain)
        el = element(domain, i)
        range = dc:dc+ndims(el)-1
        push!(ExtensionFrames, ExtensionFrame(el, element(basis, range)))
        dc += ndims(el)
    end
    tensorproduct(ExtensionFrames...)
end

extensionframe(domain, basis) = ExtensionFrame(domain, basis)

left(d::ExtensionFrame, x...) = left(domain(d))
right(d::ExtensionFrame, x...) = right(domain(d))
