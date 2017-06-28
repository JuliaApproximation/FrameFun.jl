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
import Base: ndims, length, eltype, size, push!, similar
import Base: ctranspose, inv

# Iteration and indexing
import Base: eachindex, start, next, done, getindex, unsafe_getindex,
    checkbounds

# Display
import Base: show, showcompact

import Base: promote, promote_rule, convert, promote_eltype

import Base: broadcast

## - Imports from Domains
import Domains: indomain
# - for intervals
import Domains: leftendpoint, rightendpoint
# - for mapped domains
import Domains: domain, mapping
# - for composite structures
import Domains: element, elements, nb_elements
# - for tensor products
import Domains: tensorproduct, ⊗


## - Imports from BasisFunctions
import BasisFunctions: src, dest, matrix, matrix!, apply!, apply_inplace!, numtype

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

import BasisFunctions: set_promote_eltype, superset, similar_set

import BasisFunctions: eval_set_element, eval_element, eval_expansion,
    call_set_expansion, name, in_support

import BasisFunctions: differentiate, ∂x, ∂y, ∂z, ∫∂x, ∫∂y, ∫∂z, ∫, is_compatible

import BasisFunctions: True, False, resize, promote_eltype

import BasisFunctions: show_setexpansion

import BasisFunctions: postprocess, plotgrid

import BasisFunctions: flatten

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

# from frames/extensionframe.jl
export ExtensionFrame, basis, domain, extensionframe
export Gram, DualGram, MixedGram
# from frames/sumframe.jl
export SumFrame, sumframe

# from frames/enrichedframe.jl
export EnrichedFrame

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
export Mandelbrot, JuliaSet


include("subgrid.jl")

include("domains/boundingbox.jl")
include("domains/extensions.jl")

include("frames/extensionframe.jl")
include("frames/sumframe.jl")
include("frames/enrichedframe.jl")

include("fun/basisdomains.jl")
include("fun/funs.jl")


include("approximation/fe_problem.jl")
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

include("randomgrid.jl")


end # module
