# FrameFun.jl

module FrameFun

using StaticArrays
using BasisFunctions
using RecipesBase

using Base.Cartesian

import Base: +, *, /, ==, |, &, -, \, <, <=, >, >=

import Base: intersect, union, isapprox, setdiff

import Base: length, eltype, size, push!, ctranspose, inv, similar

import Base: eachindex, start, next, done, getindex, in, unsafe_getindex,
    checkbounds

import Base: show, showcompact, call

import Base: promote, promote_rule, convert, promote_eltype

import Base: ndims, unsafe_getindex


# Imports from BasisFunctions follow
import BasisFunctions: composite_length, ⊗, tensorproduct, flatten,
    compose, elements, element, ⊕

import BasisFunctions: src, dest, matrix, matrix!, apply!, apply_inplace!, numtype

import BasisFunctions: grid, left, right, stepsize, sample

import BasisFunctions: is_basis, is_frame, is_orthogonal, is_orthonormal, is_biorthogonal,
    has_transform, has_grid, has_derivative,
    has_antiderivative, has_extension, has_grid_transform

import BasisFunctions: operator, matrix, is_diagonal, is_inplace

import BasisFunctions: coefficients, set,
    transform_operator_pre, transform_operator_post, evaluation_operator, interpolation_operator,
    differentiation_operator, antidifferentiation_operator, approximation_operator,
    extend, extension_size, extension_operator, restriction_operator,
    default_approximation_operator, has_extension, wrap_operator

import BasisFunctions: set_promote_eltype, superset, similar_set, mapping

import BasisFunctions: eval_set_element, eval_element, eval_expansion,
    call_set_expansion, name, in_support

import BasisFunctions: differentiate, ∂x, ∂y, ∂z, ∫∂x, ∫∂y, ∫∂z, ∫, is_compatible

import BasisFunctions: True, False, complexify, resize, promote_eltype

import BasisFunctions: left, right

import BasisFunctions: show_setexpansion

import BasisFunctions: postprocess, plotgrid

import BasisFunctions: Gram, DualGram, MixedGram, DiscreteGram, DiscreteDualGram, DiscreteMixedGram

import BasisFunctions: discrete_approximation_operator, continuous_approximation_operator

# about subgrids
import BasisFunctions: AbstractSubGrid, IndexSubGrid, is_subindex, supergrid,
    similar_subgrid

# from box.jl
export BBox, BBox1, BBox2, left, right
export ⊂

# from domains.jl
export Interval, Disk, Square, Cube, Ball, Cylinder, atomium, boundingbox
export ⊗, ∩, composite_length, element, elements, is_composite
export randomcircles

# from derived_domains.jl
export Characteristic
export numtype

# from funs.jl
export ExpFun, ChebyFun, Fun, SetFun, sampling_grid, domain, abserror

# from frames/extensionframe.jl
export ExtensionFrame, basis, domain, extensionframe
export Gram, DualGram, MixedGram
# from frames/sumframe.jl
export WeightedSumFrame, sumframe

# from frames/enrichedframe.jl
export EnrichedFrame

# from domains/fractal.jl
export Mandelbrot, JuliaSet

# from DiffEquation.jl
export BoundaryCondition, DiffEquation, solve

# from constructors.jl
export FunConstructor

# from space.jl
export FourierSpace, ChebyshevSpace, ⊕, add, construct
# from recipes.jl


include("box.jl")

include("domains/domains.jl")

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

include("domains/fractals.jl")


end # module
