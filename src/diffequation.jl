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
    D      :: AbstractDomain
    dRhs   :: Function
    function BoundaryCondition(S :: FunctionSet, diff :: AbstractOperator, D :: AbstractDomain, dRhs :: Function)
        new(S,diff,D,dRhs)
    end
end

BoundaryCondition(S :: FunctionSet, D::AbstractDomain) = BoundaryCondition(S,IdentityOperator(S),D,default_boundary_condition)

default_boundary_condition(x) = 0
default_boundary_condition(x,y) = 0
default_boundary_condition(x,y,z) = 0
    

immutable DiffEquation
    S     :: FunctionSet
    Diff  :: AbstractOperator
    DRhs   :: Function
    BC    :: BoundaryCondition
    D    :: AbstractDomain
    function DiffEquation(S::FunctionSet, Diff::AbstractOperator, DRhs:: Function, BC :: BoundaryCondition)
        new(S,Diff,DRhs,BC,BC.D)
    end
end


function operator(D::DiffEquation)
    problem = FE_DiscreteProblem(D.D,D.S,2)
    B = frequency_basis(problem)
    ADiff = inv(D.Diff)
    A1 = operator(problem)
    A2 = evaluation_operator(D.S,boundary(grid(src(A1)),D.D))*(D.BC.diff)*ADiff*normalization(problem)
    [A1; A2]
end

function rhs(D::DiffEquation)
    problem = FE_DiscreteProblem(D.D,D.S,2)
    op = operator(problem)
    rhs1 = sample(grid(time_basis_restricted(problem)),D.DRhs, eltype(src(op)))
    rhs2 = sample(boundary(grid(src(op)),D.D),D.BC.dRhs, eltype(src(op)))
    [rhs1; rhs2]
end

function solve(D::DiffEquation, solver=FE_ProjectionSolver; options...)
    problem = FE_DiscreteProblem(D.D,D.S,2)
    Adiff= inv(D.Diff)
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

immutable setDCZeroOperator{SRC,DEST} <: AbstractOperator{SRC,DEST}
    src :: SRC
    dest :: DEST
end

