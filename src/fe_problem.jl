# fe_problem.jl

abstract FE_Problem

# This type groups the data corresponding to a FE problem.
immutable FE_DiscreteProblem <: FE_Problem
    fbasis1         ::  AbstractBasis
    fbasis2         ::  AbstractBasis

    tbasis1         ::  AbstractBasis
    tbasis2         ::  AbstractBasis

    restricted_tbasis   ::  AbstractBasis

    f_extension     ::  AbstractOperator
    f_restriction   ::  AbstractOperator

    t_extension     ::  AbstractOperator
    t_restriction   ::  AbstractOperator

    transform1      ::  AbstractOperator
    itransform1     ::  AbstractOperator

    transform2      ::  AbstractOperator
    itransform2     ::  AbstractOperator

    op              ::  AbstractOperator
    opt             ::  AbstractOperator

    function FE_DiscreteProblem(fbasis1, fbasis2, tbasis1, tbasis2, restricted_tbasis, 
        f_extension, f_restriction, t_extension, t_restriction, 
        transform1, itransform1, transform2, itransform2)

        op  = t_restriction * itransform2 * f_extension
        opt = f_restriction * transform2 * t_extension

        new(fbasis1, fbasis2, tbasis1, tbasis2, restricted_tbasis, 
            f_extension, f_restriction, t_extension, t_restriction, 
            transform1, itransform1, transform2, itransform2, 
            op, opt)
    end

end


numtype(p::FE_DiscreteProblem) = numtype(p.fbasis1)

eltype(p::FE_DiscreteProblem) = eltype(operator(p))

dim(p::FE_DiscreteProblem) = dim(p.fbasis1)

size(p::FE_DiscreteProblem) = size(operator(p))

size(p::FE_DiscreteProblem, j) = size(operator(p), j)

frequency_basis(p::FE_DiscreteProblem) = p.fbasis1

frequency_basis_ext(p::FE_DiscreteProblem) = p.fbasis2

time_basis(p::FE_DiscreteProblem) = p.tbasis1

time_basis_ext(p::FE_DiscreteProblem) = p.tbasis2

restricted_time_basis(p::FE_DiscreteProblem) = p.restricted_tbasis

operator(p::FE_DiscreteProblem) = p.op

operator_transpose(p::FE_DiscreteProblem) = p.opt


function rhs(p::FE_DiscreteProblem, f::Function, elt = eltype(op))
    grid1 = grid(restricted_time_basis(p))
    M = length(grid1)
    b = Array(elt, M)
    rhs!(p, b, f)
    b
end

function rhs!(p::FE_DiscreteProblem, b::Vector, f::Function)
    grid1 = grid(restricted_time_basis(p))
    M = length(grid1)

    @assert length(b) == M

    rhs!(grid1, b, f)
end

function rhs!(grid::AbstractGrid, b::Vector, f::Function)
    l = 0
    for i in eachindex(grid)
        l = l+1
        b[l] = f(grid[i])
    end
end

function rhs!(grid::MaskedGrid, b::Vector, f::Function)
    l = 0
    for i in eachindex_mask(grid)
        if in(i, grid)
            l = l+1
            b[l] = f(grid[i])
        end
    end
end



