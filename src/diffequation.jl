# diffequation.jl

""" 
A DiffEquation describes a differential equation, with or without boundary conditions.
Parameters:
- Fun is the SetFun that will describe the result

When the equation is solved the equations:
-Diff(Fun) = DRhs on the domain
-diff(Fun) = dRhs on boundary
will hold.
"""

immutable BoundaryCondition
    S      :: FunctionSet
    diff   :: AbstractOperator
    DG      :: AbstractGrid
    dRhs   :: Function
    function BoundaryCondition(S :: FunctionSet, diff :: AbstractOperator, DG :: AbstractGrid, dRhs :: Function)
        new(S,diff,DG,dRhs)
    end
end

BoundaryCondition(S :: FunctionSet, D::AbstractDomain) = BoundaryCondition(S,IdentityOperator(S),boundary(grid(S),D),default_boundary_condition)
BoundaryCondition(S :: FunctionSet, diff::AbstractOperator, D::AbstractDomain) = BoundaryCondition(S,diff,boundary(grid(S),D),default_boundary_condition)
BoundaryCondition(S :: FunctionSet, diff::AbstractOperator, D::AbstractDomain, dRhs::Function) = BoundaryCondition(S,diff,boundary(grid(S),D),dRhs)
default_boundary_condition(x) = 0
default_boundary_condition(x,y) = 0
default_boundary_condition(x,y,z) = 0
    

immutable DiffEquation
    S     :: FunctionSet
    D     :: AbstractDomain
    Diff  :: AbstractOperator
    DRhs   :: Function
    BCs    :: Tuple
    function DiffEquation(S::FunctionSet, D::AbstractDomain,Diff::AbstractOperator, DRhs:: Function, BCs::Tuple)
        new(S,D,Diff,DRhs,BCs)
    end
end

DiffEquation(S::FunctionSet, D::AbstractDomain, Diff::AbstractOperator, DRhs::Function, BC::BoundaryCondition) = DiffEquation(S,D,Diff,DRhs,(BC,))

function operator(D::DiffEquation)
    problem = FE_DiscreteProblem(D.D,D.S,2)
    B = frequency_basis(problem)
    ADiff = inv_diagonal(D.Diff)
    ops = Array{AbstractOperator}(length(D.BCs)+1,1)
    ops[1] = operator(problem)*D.Diff*ADiff
    for i = 1:length(D.BCs)
        Ac = evaluation_operator(D.S,D.BCs[i].DG)*(D.BCs[i].diff)*ADiff*normalization(problem)
        ops[i+1]=Ac
    end
    BlockOperator(ops)
end

function rhs(D::DiffEquation)
    problem = FE_DiscreteProblem(D.D,D.S,2)
    op = operator(problem)
    rhs = Array(Array{eltype(src(op)),1},0)
    push!(rhs,sample(grid(time_basis_restricted(problem)),D.DRhs, eltype(src(op))))
    for BC in D.BCs
        push!(rhs,sample(BC.DG,BC.dRhs, eltype(src(op))))
    end
    MultiArray(rhs)
end

function solve(D::DiffEquation, solver=FE_ProjectionSolver; options...)
    Adiff= inv_diagonal(D.Diff)
    b = rhs(D)
    DEproblem = problem(D)
    A = solver(DEproblem; options...)
    coef  = A * b
    SetFun(D.D, dest(A), Adiff*coef)
end

function problem(D::DiffEquation)
    problem = FE_DiscreteProblem(D.D,D.S,2)
    op = operator(D)
    opt = ctranspose(op)
    fb = frequency_basis(problem)
    fbe = frequency_basis_ext(problem)
    tb = MultiSet([time_basis(problem); elements(dest(op))[2:end]])
    tbe = MultiSet([time_basis_ext(problem); elements(dest(op))[2:end]])
    tbr = MultiSet([time_basis_restricted(problem); elements(dest(op))[2:end]])
    fe = f_extension(problem)
    fr = f_restriction(problem)
    te = t_extension(problem)⊕IdentityOperator(element(dest(op),2:length(elements(dest(op)))))
    tr = t_restriction(problem)⊗IdentityOperator(element(dest(op),2:length(elements(dest(op)))))
    DEproblem = FE_DiscreteProblem(domain(problem),op, opt, fb,fbe,tb,tbe,tbr,fe,fr,te,tr, transform1(problem), itransform1(problem), transform2(problem), itransform2(problem),normalization(problem))
end



