# FrameFuns.jl

module FrameFuns

using FixedSizeArrays
using BasisFunctions
#using PyPlot
using RecipesBase
using Compat

using Base.Cartesian

import Base: +, *, /, ==, |, &, -, \, <, <=, >, >=

import Base: intersect, union, isapprox, setdiff

import Base: length, eltype, size, push!, ctranspose, similar

import Base: eachindex, start, next, done, getindex, in

import Base: show, showcompact, call

import Base: promote, promote_rule, convert, promote_eltype

import Base: ndims

# import PyPlot: plot

import BasisFunctions: composite_length, ⊗, tensorproduct, flatten,
    compose, elements, element, ⊕

import BasisFunctions: src, dest, matrix, matrix!, apply!, apply_inplace!, numtype

import BasisFunctions: grid, left, right, stepsize, sample

import BasisFunctions: operator, coefficients, set, is_basis, is_frame, is_diagonal, is_inplace,
    transform_pre_operator, transform_post_operator, evaluation_operator, interpolation_operator,
    differentiation_operator, antidifferentiation_operator, approximation_operator,
    extend, extension_size, extension_operator, restriction_operator,
    default_approximation_operator, has_extension, wrap_operator

import BasisFunctions: call_set, call_set!, call_expansion_with_set,
call_expansion_with_set!, call_expansion, call_expansion!, call_element, name

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

# from derived_domains.jl
export Characteristic
export numtype

# from funs.jl
export ExpFun, ChebyFun, Fun, FrameFun, sampling_grid, domain, abserror

# from domainframe.jl
export DomainFrame, basis, call_set, call_set!

# from fractal.jl
export Mandelbrot, JuliaSet

# from DiffEquation.jl
export BoundaryCondition, DiffEquation, solve

# from constructors.jl
export FunConstructor

# from space.jl
export FourierSpace, ChebyshevSpace, ⊕, add, construct

# from plots.jl
#export plot, plot_error, plot_samples, plot_domain, plot_image
# from recipes.jl

# We support both vectors (AbstractVector) and FixedSizeArray's (Vec)
typealias AnyVector Union{AbstractVector,Vec}


include("box.jl")

include("domains.jl")

include("fractals.jl")

include("subgrid.jl")

include("domainframe.jl")

include("funs.jl")

include("fourierdomains.jl")

include("fe_problem.jl")

include("fe_solvers.jl")

include("fastsolver.jl")

include("smoothsolver.jl")

include("fe_approx.jl")

include("recipes.jl")

include("diffequation.jl")

include("space.jl")

include("constructors.jl")

end # module
