# concat_solver.jl

"""
Solver to approximate functions of the form f(x) = p(x) + w(x) q(x), where w is a given
weight function, and p and q are well approximated by set1 and set2 respectively.
"""
immutable ConcatSolver{ELT} <: AbstractOperator{ELT}
    src     ::  FunctionSet
    dest    ::  FunctionSet
    set1    ::  FunctionSet
    set2    ::  FunctionSet
    fun

    A1
    A2
    D
    P
    C
    A1prime
    A2prime
    normalization1
    normalization2
    lr
end


# function ConcatSolver(concatset, set1, set2, weightfunction)
#     ELT = eltype(concatset)
#     extended_set = resize(set1, 2*(length(set1)+length(set2)))
#     ls_grid = grid(extended_set)
#     dgs = DiscreteGridSpace(extended_set)
#     T = transform_operator(extended_set, dgs)
#     N = transform_normalization_operator(extended_set)
#     E1 = extension_operator(set1, extended_set)
#     E2 = extension_operator(set2, extended_set)
#     R1 = restriction_operator(extended_set, set1)
#     R2 = restriction_operator(extended_set, set2)
#     D = matrix(DiagonalOperator(dgs, ELT[weightfunction(x) for x in ls_grid]))
#     A1 = matrix(T*E1)
#     A2 = matrix(T*E2)
#     P = eye(ELT,length(extended_set)) - A1*A1'
#     C = P*D*A2
#     normalization1 = matrix(R1*N*E1)
#     normalization2 = matrix(R2*N*E2)
#     ConcatSolver{ELT}(dgs, concatset, set1, set2, weightfunction,
#         A1, A2, D, P, C, A1', A2', normalization1, normalization2)
# end


function ConcatSolver(concatset, set1, set2, weightfunction)
    ELT = eltype(concatset)
    extended_set = resize(set1, 2*(length(set1)+length(set2)))
    ls_grid = grid(extended_set)
    dgs = DiscreteGridSpace(extended_set)
    T = transform_operator(extended_set, dgs)
    N = transform_normalization_operator(extended_set)
    E1 = extension_operator(set1, extended_set)
    E2 = extension_operator(set2, extended_set)
    R1 = restriction_operator(extended_set, set1)
    R2 = restriction_operator(extended_set, set2)
    D = DiagonalOperator(dgs, ELT[weightfunction(x) for x in ls_grid])
    A1 = T*E1
    A2 = T*E2
    P = IdentityOperator(dgs) - A1*A1'
    C = P*D*A2
    lr = lowrank_approximation(C)
    normalization1 = matrix(R1*N*E1)
    normalization2 = matrix(R2*N*E2)
    ConcatSolver{ELT}(dgs, concatset, set1, set2, weightfunction,
        A1, A2, D, P, C, A1', A2', normalization1, normalization2, lr)
end


function apply!(op::ConcatSolver, coef_dest, coef_src)
    B = coef_src
    PB = op.P * B
    x2 = lowranksolve(op.lr, PB)
    B2 = B - op.D * (op.A2 * x2)
    x1 = op.A1prime * B2
    x1b = op.normalization1*x1
    x2b = op.normalization2*x2
    coef_dest[:] = [x1b; x2b]
end

# function apply!(op::ConcatSolver, coef_dest, coef_src)
#     B = coef_src
#     PB = op.P * B
#     x2 = op.C \ PB
#     B2 = B - op.D * op.A2 * x2
#     x1 = op.A1prime * B2
#     x1b = op.normalization1*x1
#     x2b = op.normalization2*x2
#     coef_dest[:] = [x1b; x2b]
# end
#

function approximation_operator(set::MultiSet; options...)
    @assert composite_length(set) == 2
    approximation_operator_concat(set, element(set, 1), element(set, 2); options...)
end

approximation_operator_concat(set::MultiSet, set1::FunctionSet, set2::FunctionSet) =
    default_approximation_operator(set)

function approximation_operator_concat(concatset::MultiSet, set1::FunctionSet, set2::AugmentedSet)
    if has_transform(set1) && has_transform(set(set2))
        ConcatSolver(concatset, set1, set(set2), fun(set2))
    else
        default_approximation_operator(concatset)
    end
end
