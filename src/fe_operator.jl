# fe_operator.jl

abstract FE_Operator{SRC,DEST} <: AbstractOperator{SRC,DEST}

numtype(op::FE_Operator) = numtype(fe_problem(op))


frequency_basis(op::FE_Operator) = frequency_basis(fe_problem(op))


immutable FE_DiscreteOperator{SRC,DEST} <: FE_Operator{SRC,DEST}
    problem ::  FE_DiscreteProblem
end

FE_DiscreteOperator(problem::FE_DiscreteProblem) = FE_DiscreteOperator{typeof(frequency_basis(problem)),typeof(restricted_time_basis(problem))}(problem)

src(op::FE_DiscreteOperator) = frequency_basis(op.problem)

dest(op::FE_DiscreteOperator) = restricted_time_basis(op.problem)

fe_problem(op::FE_DiscreteOperator) = op.problem


function apply!(op::FE_DiscreteOperator, coef_dest, coef_src)
    p = op.problem
    apply!(p.f_extension, p.scratch1, coef_src)
    apply!(p.itransform2, p.scratch2, p.scratch1)
    apply!(p.t_restriction, coef_dest, p.scratch2)
end


function fe_rhs(op::FE_DiscreteOperator, f::Function, elt = eltype(op))
    grid = natural_grid(restricted_time_basis(fe_problem(op)))
    M = length(grid)
    b = Array(elt, M)
    fe_rhs!(op, b, f)
    b
end

function fe_rhs!(op::FE_DiscreteOperator, b::Vector, f::Function)
    grid = natural_grid(restricted_time_basis(fe_problem(op)))
    M = length(grid)

    @assert length(b) == M

    fe_rhs!(grid, b, f)
end

function fe_rhs!(grid::AbstractGrid, b::Vector, f::Function)
    l = 0
    for i in eachindex(grid)
        l = l+1
        b[l] = f(grid[i])
    end
end

function fe_rhs!(grid::MaskedGrid, b::Vector, f::Function)
    l = 0
    for i in eachindex_mask(grid)
        if in(i, grid)
            l = l+1
            b[l] = f(grid[i])
        end
    end
end


