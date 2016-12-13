# FrameFun.jl

module FrameFun

using StaticArrays
using BasisFunctions
using RecipesBase

using Base.Cartesian

import Base: +, *, /, ==, |, &, -, \, <, <=, >, >=

import Base: intersect, union, isapprox, setdiff

import Base: length, eltype, size, push!, ctranspose, inv, similar

import Base: eachindex, start, next, done, getindex, in

import Base: show, showcompact, call

import Base: promote, promote_rule, convert, promote_eltype

import Base: ndims


# Imports from BasisFunctions follow
import BasisFunctions: composite_length, ⊗, tensorproduct, flatten,
    compose, elements, element, ⊕

import BasisFunctions: src, dest, matrix, matrix!, apply!, apply_inplace!, numtype

import BasisFunctions: grid, left, right, stepsize, sample

import BasisFunctions: operator, coefficients, set, is_basis, is_frame, is_diagonal, is_inplace,
    transform_operator_pre, transform_operator_post, evaluation_operator, interpolation_operator,
    differentiation_operator, antidifferentiation_operator, approximation_operator,
    extend, extension_size, extension_operator, restriction_operator,
    default_approximation_operator, has_extension, wrap_operator

import BasisFunctions: eval_set_element, eval_element, eval_expansion,
    call_set_expansion, name

import BasisFunctions: differentiate, ∂x, ∂y, ∂z, ∫∂x, ∫∂y, ∫∂z, ∫, is_compatible

import BasisFunctions: True, False, complexify, resize, promote_eltype

import BasisFunctions: left, right

import BasisFunctions: show_setexpansion

import BasisFunctions: postprocess, plotgrid

# about subgrids
import BasisFunctions: AbstractSubGrid, IndexSubGrid, is_subindex, supergrid,
    similar_subgrid

# from box.jl
export BBox, BBox1, BBox2, left, right
export ⊂

# from domains.jl
export Interval, Disk, Square, Cube, Ball, Cylinder, atomium, boundingbox
export ⊗, ∩, composite_length, element, elements

# from derived_domains.jl
export Characteristic
export numtype

# from funs.jl
export ExpFun, ChebyFun, Fun, SetFun, sampling_grid, domain, abserror

# from domainframe.jl
export DomainFrame, basis

# from fractal.jl
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

include("frames/domainframe.jl")

include("fun/funs.jl")

include("approximation/fe_problem.jl")

include("approximation/fe_solvers.jl")

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
