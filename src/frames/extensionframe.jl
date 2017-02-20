# extensionframe.jl

"""
An ExtensionFrame is the restriction of a basis to a subset of its domain. This
results in a frame that implicitly represents extensions of functions on the
smaller set to the larger set.
"""
immutable ExtensionFrame{N,T} <: DerivedSet{N,T}
    domain      ::  AbstractDomain{N}
    basis       ::  FunctionSet{N,T}

    function ExtensionFrame(domain::AbstractDomain, basis::FunctionSet)
        @assert is_basis(basis)
        new(domain, basis)
    end
end

ExtensionFrame{N,T}(domain::AbstractDomain{N}, basis::FunctionSet{N,T}) =
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
function extensionframe(domain::TensorProductDomain, basis::TensorProductSet)
    ExtensionFrames = FunctionSet[]
    dc = 1
    for i = 1:composite_length(domain)
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

# function Gram(f::ExtensionFrame; options...)
#   A = zeros(eltype(f),length(f),length(f))
#   grammatrix!(A,f; options...)
#   MatrixOperator(A)
# end

#TODO DualGram executed is of functionset.jl but the is_biorthogonal trait is not correct

MixedGram(f::ExtensionFrame; options...) = DualGram(basis(f))*Gram(f; options...)

function grammatrix!(result, frame::ExtensionFrame; options...)
  b = basis(frame)
  for i in 1:size(result,1)
    for j in i:size(result,2)
      I = innerproduct(frame, x->b[i](x), j; options...)
      result[i,j] = I
      if i!= j
        result[j,i] = I
      end
    end
  end
  result
end

import BasisFunctions: innerproduct
innerproduct(frame::ExtensionFrame, f::Function, idx::Int; options...) =
    innerproduct(basis(frame), domain(frame), f, idx; options...)

# TODO now we assume that domainunion contains sections that do not overlap
function innerproduct(b::FunctionSet, d::DomainUnion, f::Function, idx::Int; options...)
  d1 = firstelement(domain)
  d2 = firstelement(domain)
  innerproduct(b, d1, f, idx; options...) + innerproduct(b, d2, f, idx; options...)
end

innerproduct(b::FunctionSet1d, d::Interval, f::Function, idx::Int; options...) =
    innerproduct(b, f, idx, left(d), right(d); options...)

import BasisFunctions: continuous_approximation_operator
continuous_approximation_operator(frame::ExtensionFrame; solver = ContinuousDirectSolver, options...) = solver(frame; options...)
