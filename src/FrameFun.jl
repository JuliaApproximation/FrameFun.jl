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


## - Imports from Domains
import Domains: left, right
import Domains: element, elements, nb_elements
import Domains: tensorproduct, flatten, ⊗
import Domains: boundingbox


## - Imports from BasisFunctions
import BasisFunctions: src, dest, matrix, matrix!, apply!, apply_inplace!, numtype

import BasisFunctions: grid, left, right, stepsize, sample

import BasisFunctions: is_basis, is_frame, has_transform, has_grid, has_derivative,
    has_antiderivative, has_extension, has_grid_transform

import BasisFunctions: operator, matrix, is_diagonal, is_inplace, ⊕

import BasisFunctions: coefficients, set,
    transform_operator_pre, transform_operator_post, evaluation_operator, interpolation_operator,
    differentiation_operator, antidifferentiation_operator, approximation_operator,
    extend, extension_size, extension_operator, restriction_operator,
    default_approximation_operator, has_extension, wrap_operator

import BasisFunctions: set_promote_eltype, superset, similar_set, mapping

import BasisFunctions: eval_set_element, eval_element, eval_expansion,
    call_set_expansion, name, in_support

import BasisFunctions: differentiate, ∂x, ∂y, ∂z, ∫∂x, ∫∂y, ∫∂z, ∫, is_compatible

import BasisFunctions: True, False, resize, promote_eltype

import BasisFunctions: show_setexpansion

import BasisFunctions: postprocess, plotgrid

# about subgrids
import BasisFunctions: AbstractSubGrid, IndexSubGrid, is_subindex, supergrid,
    similar_subgrid


###############################
## Exhaustive list of exports
###############################

# from funs.jl
export ExpFun, ChebyFun, Fun, SetFun, sampling_grid, domain, abserror

# from subgrid.jl
export MaskedGrid

# from frames/extensionframe.jl
export ExtensionFrame, basis, domain

# from frames/sumframe.jl
export SumFrame, sumframe

# from frames/enrichedframe.jl
export EnrichedFrame

# from DiffEquation.jl
export BoundaryCondition, DiffEquation, solve

# from constructors.jl
export FunConstructor

# from space.jl
export FourierSpace, ChebyshevSpace, ⊕, add, construct
# from recipes.jl

# from randomgrid.jl
export randomgrid, randompoint


include("subgrid.jl")

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

include("recipes.jl")

include("diffequation.jl")

include("approximation/space.jl")
include("approximation/constructors.jl")

include("domains/fourierdomains.jl")

include("randomgrid.jl")


end # module
