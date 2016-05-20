# FrameFuns.jl

module FrameFuns

using FixedSizeArrays
using BasisFunctions
using Debug
#using ApproxFun

using Base.Cartesian

import Base: +, *, /, ==, |, &, -, \

import Base: intersect, union, isapprox, setdiff

import Base: length, eltype, size, push!, ctranspose, similar

import Base: eachindex, start, next, done, getindex, in

import Base: show, showcompact, call, convert

import Base: ndims

import BasisFunctions: composite_length, ⊗, tensorproduct, flatten,
    compose, elements, element

import BasisFunctions: src, dest, matrix, matrix!, apply!, numtype

import BasisFunctions: index_dim, grid, left, right, stepsize, sample

import BasisFunctions: operator, coefficients, set, is_basis, is_frame,
    transform_normalization_operator, evaluation_operator, interpolation_operator,
    differentiation_operator, antidifferentiation_operator, approximation_operator,
    extend, extension_size, default_approximation_operator

import BasisFunctions: call_set, call_set!, call_expansion, call_expansion!, call_element, name

import BasisFunctions: differentiate, ∂x, ∂y, ∂z, ∫∂x, ∫∂y, ∫∂z, ∫

import BasisFunctions: True, False, complexify, resize, promote_eltype

import BasisFunctions: left, right

import BasisFunctions: show_setexpansion

# from box.jl
export BBox, BBox1, BBox2, left, right
export ⊂

# from domains.jl
export Interval, Disk, Square, Cube, Ball, Cylinder, atomium, boundingbox
export ⊗, ∩, composite_length, element, elements

export numtype

# from funs.jl
export ExpFun, ChebyFun, Fun, FrameFun, sampling_grid, domain

# from domainframe.jl
export DomainFrame, basis, call_set, call_set!

# from fractal.jl
export Mandelbrot, JuliaSet

# from plots.jl
export plot, plot_expansion, plot_error, plot_samples, plot_domain, plot_image

# We support both vectors (AbstractVector) and FixedSizeArray's (Vec)
typealias AnyVector Union{AbstractVector,Vec}


include("box.jl")

include("domains.jl")

include("fractals.jl")

include("subgrid.jl")

include("domainframe.jl")

include("funs.jl")

include("lowrank.jl")

include("fe_problem.jl")

include("fe_solvers.jl")

include("fastsolver.jl")

include("smoothsolver.jl")

include("fe_approx.jl")

include("concat_solver.jl")

include("explicit_extension.jl")

# TODO: try out Plots.jl
include("plots.jl")

end # module
