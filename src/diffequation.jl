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
    A1 = operator(problem)
    A2 = evaluation_operator(D.S,boundary(grid(src(A1)),D.D))*(D.BC.diff)
    [A1; A2]
end

function Rhs(D::DiffEquation)
    problem = FE_DiscreteProblem(D.D,D.S,2)
    op = operator(problem)
    rhs1 = sample(grid(time_basis_restricted(problem)),D.DRhs, eltype(src(op)))
    rhs2 = sample(boundary(grid(src(op)),D.D),D.BC.dRhs, eltype(src(op)))
    [rhs1; rhs2]
end
        
