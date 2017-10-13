# FrameFun.jl

module FrameFun

using Base.Cartesian
using StaticArrays
using RecipesBase

using Domains
using BasisFunctions

###############################
## Exhaustive list of imports
###############################

## - Imports from Base
import Base: +, *, /, ^, ==, |, &, -, \, <, <=, >, >=

import Base: intersect, union, isapprox, setdiff, in

# Arrays
import Base: length, eltype, size, push!, similar
import Base: ctranspose, inv

# Iteration and indexing
import Base: eachindex, start, next, done, getindex, unsafe_getindex,
    checkbounds

# Display
import Base: show, showcompact

import Base: promote, promote_rule, convert, promote_eltype

import Base: broadcast

## - Imports from Domains
import Domains: indomain, dimension
# - for intervals
import Domains: leftendpoint, rightendpoint
# - for mapped domains
import Domains: domain
# - for composite structures
import Domains: element, elements, nb_elements
# - for cartesian products
import Domains: cartesianproduct, ×


## - Imports from BasisFunctions
import BasisFunctions: src, dest, matrix, matrix!, apply!, apply_inplace!, dimension, rangetype, domaintype

import BasisFunctions: tensorproduct, ⊗

import BasisFunctions: grid, left, right, stepsize, sample

import BasisFunctions: is_basis, is_frame, is_orthogonal, is_orthonormal, is_biorthogonal,
    has_transform, has_grid, has_derivative,
    has_antiderivative, has_extension, has_grid_transform

import BasisFunctions: operator, matrix, is_diagonal, is_inplace, ⊕

import BasisFunctions: coefficients, set,
    transform_operator_pre, transform_operator_post, evaluation_operator, interpolation_operator,
    differentiation_operator, antidifferentiation_operator, approximation_operator,
    extend, extension_size, extension_operator, restriction_operator,
    default_approximation_operator, has_extension, wrap_operator

import BasisFunctions: superset, similar_set, promote_domaintype, promote_domainsubtype

import BasisFunctions: eval_set_element, eval_element, eval_expansion,
    call_set_expansion, name, in_support

import BasisFunctions: differentiate, ∂x, ∂y, ∂z, ∫∂x, ∫∂y, ∫∂z, ∫, is_compatible

import BasisFunctions: True, False, resize, promote_eltype

import BasisFunctions: show_setexpansion

import BasisFunctions: postprocess, plotgrid

import BasisFunctions: flatten

import BasisFunctions: Span, Span1d, Span2d

import BasisFunctions: span, coefficient_type, coeftype, similar_span

# about subgrids
import BasisFunctions: AbstractSubGrid, IndexSubGrid, is_subindex, supergrid,
    similar_subgrid

import BasisFunctions: Gram, DualGram, MixedGram, DiscreteGram, DiscreteDualGram, DiscreteMixedGram

import BasisFunctions: discrete_approximation_operator, continuous_approximation_operator

###############################
## Exhaustive list of exports
###############################
# from funs.jl
export ExpFun, ChebyFun, Fun, SetFun, sampling_grid, domain, abserror

# from subgrid.jl
export MaskedGrid

# from domains/boundingbox.jl
export BoundingBox, BBox, BBox1, BBox2, BBox3, BBox4
export boundingbox

# from domains/extensions.jl
export dist, normal

# from frames/extensionframe.jl
export ExtensionFrame, basis, domain, extensionframe
export Gram, DualGram, MixedGram

# from DiffEquation.jl
export DirichletBC, NeumannBC, DiffEquation, solve

# from constructors.jl
export FunConstructor

# from space.jl
export FourierSpace, ChebyshevSpace, ⊕, add, construct
# from recipes.jl

# from randomgrid.jl
export randomgrid, randompoint

# from domains/fractals.jl
export mandelbrot, juliaset

# from domains/atomium.jl
export atomium

# from domains/polardomain.jl
export polardomain

# from domains/characteristic.jl
export characteristic

include("subgrid.jl")

#include("domains/boundingbox.jl")
include("domains/extensions.jl")

include("frames/extensionframe.jl")

include("fun/basisdomains.jl")
include("fun/funs.jl")


include("approximation/fe_solvers.jl")
include("approximation/lowranksolver.jl")
include("approximation/fastsolver.jl")
include("approximation/smoothsolver.jl")
include("approximation/fe_approx.jl")

include("approximation/oversampling.jl")

include("recipes.jl")

include("diffequation.jl")

include("approximation/space.jl")
include("approximation/constructors.jl")

include("domains/fourierdomains.jl")
include("domains/fractals.jl")
include("domains/atomium.jl")
include("domains/characteristic.jl")
include("domains/polardomain.jl")

include("randomgrid.jl")
include("oversampledgrid.jl")


end # module
