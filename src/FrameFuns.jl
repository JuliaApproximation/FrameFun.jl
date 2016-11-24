# FrameFuns.jl

module FrameFuns

using StaticArrays
using BasisFunctions
using RecipesBase

using Base.Cartesian

import Base: +, *, /, ==, |, &, -, \

import Base: intersect, union, isapprox, setdiff

import Base: length, eltype, size, push!, ctranspose, inv, similar

import Base: eachindex, start, next, done, getindex, in

import Base: show, showcompact, call, convert

import Base: ndims


# Imports from BasisFunctions follow
import BasisFunctions: composite_length, ⊗, tensorproduct, flatten,
    compose, elements, element

import BasisFunctions: src, dest, matrix, matrix!, apply!, apply_inplace!, numtype

import BasisFunctions: grid, left, right, stepsize, sample

import BasisFunctions: operator, coefficients, set, is_basis, is_frame, is_diagonal, is_inplace,
    transform_operator_pre, transform_operator_post, evaluation_operator, interpolation_operator,
    differentiation_operator, antidifferentiation_operator, approximation_operator,
    extend, extension_size, extension_operator, restriction_operator,
    default_approximation_operator, has_extension

import BasisFunctions: eval_set_element, eval_element, eval_expansion,
    call_set_expansion, name

import BasisFunctions: differentiate, ∂x, ∂y, ∂z, ∫∂x, ∫∂y, ∫∂z, ∫, is_compatible

import BasisFunctions: True, False, complexify, resize, promote_eltype

import BasisFunctions: left, right

import BasisFunctions: show_setexpansion

import BasisFunctions: postprocess, plotgrid

# from box.jl
export BBox, BBox1, BBox2, left, right
export ⊂

# from domains.jl
export Interval, Disk, Square, Cube, Ball, Cylinder, atomium, boundingbox
export ⊗, ∩, composite_length, element, elements

export numtype

# from funs.jl
export ExpFun, ChebyFun, Fun, FrameFun, sampling_grid, domain, abserror

# from domainframe.jl
export DomainFrame, basis

# from fractal.jl
export Mandelbrot, JuliaSet

# from DiffEquation.jl
export BoundaryCondition, DiffEquation, solve

# from recipes.jl


include("box.jl")

include("domains.jl")

include("fractals.jl")

include("subgrid.jl")

include("domainframe.jl")

include("funs.jl")

include("fe_problem.jl")

include("fe_solvers.jl")

include("fastsolver.jl")

include("smoothsolver.jl")

include("fe_approx.jl")

include("recipes.jl")

include("diffequation.jl")

end # module
