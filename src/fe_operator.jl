# fe_operator.jl

abstract FE_Operator{SRC,DEST} <: AbstractOperator{SRC,DEST}

numtype(op::FE_Operator) = numtype(problem(op))


frequency_basis(op::FE_Operator) = frequency_basis(problem(op))


immutable FE_DiscreteOperator{SRC,DEST} <: FE_Operator{SRC,DEST}
    problem ::  FE_DiscreteProblem

    scratch1
    scratch2

    function FE_DiscreteOperator(problem)
        scratch1 = Array(eltype(problem), size(frequency_basis_ext(problem)))
        scratch2 = Array(eltype(problem), size(frequency_basis_ext(problem)))
        new(problem, scratch1, scratch2)
    end
end

FE_DiscreteOperator(problem::FE_DiscreteProblem) = FE_DiscreteOperator{typeof(frequency_basis(problem)),typeof(restricted_time_basis(problem))}(problem)

src(op::FE_DiscreteOperator) = frequency_basis(op.problem)

dest(op::FE_DiscreteOperator) = restricted_time_basis(op.problem)

problem(op::FE_DiscreteOperator) = op.problem


function apply!(op::FE_DiscreteOperator, coef_dest, coef_src)
    p = problem(op)
    apply!(p.f_extension, op.scratch1, coef_src)
    apply!(p.itransform2, op.scratch2, op.scratch1)
    apply!(p.t_restriction, coef_dest, op.scratch2)
end


function rhs(op::FE_DiscreteOperator, f::Function, elt = eltype(op))
    grid = natural_grid(restricted_time_basis(problem(op)))
    M = length(grid)
    b = Array(elt, M)
    rhs!(op, b, f)
    b
end

function rhs!(op::FE_DiscreteOperator, b::Vector, f::Function)
    grid = natural_grid(restricted_time_basis(problem(op)))
    M = length(grid)

    @assert length(b) == M

    rhs!(grid, b, f)
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


