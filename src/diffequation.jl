# diffequation.jl

""" 
A DiffEquation describes a differential equation, with or without boundary conditions.
Parameters:
- Fun is the FrameFun that will describe the result

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
    A = operator(problem)*D.Diff*ADiff
    Ac = evaluation_operator(D.S,D.BCs[1].DG)*(D.BCs[1].diff)*ADiff*normalization(problem)
    for i = 2:length(D.BCs)
        Abc = evaluation_operator(D.S,D.BCs[i].DG)*(D.BCs[i].diff)*ADiff*normalization(problem)
        Ac = [Ac; Abc]
    end
    A = [A; Ac]
end

function rhs(D::DiffEquation)
    problem = FE_DiscreteProblem(D.D,D.S,2)
    op = operator(problem)
    rhs1 = sample(grid(time_basis_restricted(problem)),D.DRhs, eltype(src(op)))
    for BC in D.BCs
        rhs2 = sample(BC.DG,BC.dRhs, eltype(src(op)))
        rhs1=[rhs1; rhs2]
    end
    rhs1
end

function solve(D::DiffEquation, solver=FE_ProjectionSolver; options...)
    problem = FE_DiscreteProblem(D.D,D.S,2)
    Adiff= inv_diagonal(D.Diff)
    op = operator(D)
    b = rhs(D)
    opt = ctranspose(op)
    DEproblem = FE_DiscreteProblem(domain(problem),op, opt, frequency_basis(problem), frequency_basis_ext(problem), time_basis(problem)⊕dest(op.op2), time_basis_ext(problem)⊕dest(op.op2),time_basis_restricted(problem)⊕dest(op.op2), f_extension(problem), f_restriction(problem), t_extension(problem)⊕IdentityOperator(dest(op.op2)), t_restriction(problem)⊗IdentityOperator(dest(op.op2)), transform1(problem), itransform1(problem), transform2(problem), itransform2(problem),normalization(problem))
    A = solver(DEproblem; options...)
    coef  = A * b
    FrameFun(D.D, dest(A), Adiff*coef, A)
end

function problem(problem::FE_DiscreteProblem, op::AbstractOperator, opt::AbstractOperator)
    fb = frequency_basis(problem)
    fbe = frequency_basis_ext(problem)
    tb = time_basis(problem)⊕dest(op.op2)
    tbe = time_basis_ext(problem)⊕dest(op.op2)
    tbr = time_basis_restricted(problem)⊕dest(op.op2)
    fe = f_extension(problem)
    fr = f_restriction(problem)
    te = t_extension(problem)⊕IdentityOperator(dest(op.op2))
    tr = t_restriction(problem)⊗IdentityOperator(dest(op.op2))
    DEproblem = FE_DiscreteProblem(domain(problem),op, opt, fb,fbe,tb,tbe,tbr,fe,fr,te,tr, transform1(problem), itransform1(problem), transform2(problem), itransform2(problem),normalization(problem))
end



