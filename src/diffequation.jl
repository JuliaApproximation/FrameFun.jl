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
    sampling_factor
    function DiffEquation(S::FunctionSet, D::AbstractDomain,Diff::AbstractOperator, DRhs:: Function, BCs::Tuple, sampling_factor=2)
        new(S,D,Diff,DRhs,BCs, sampling_factor)
    end
end

DiffEquation(S::FunctionSet, D::AbstractDomain, Diff::AbstractOperator, DRhs::Function, BC::BoundaryCondition, sampling_factor=2) = DiffEquation(S,D,Diff,DRhs,(BC,), sampling_factor)

function operator(D::DiffEquation; incboundary=false, options...)
    #problem = FE_DiscreteProblem(D.D,D.S,2)
    B = D.S
    ADiff = inv_diagonal(D.Diff)
    ops = incboundary ? Array{AbstractOperator}(length(D.BCs)+2,1) : Array{AbstractOperator}(length(D.BCs)+1,1)
    G, lB = oversampled_grid(D.D,D.S,D.sampling_factor)
    
    op = grid_evaluation_operator(D.S,DiscreteGridSpace(G,eltype(D.S)),G)
    
    cnt=1
    ops[cnt] = op*D.Diff*ADiff
    if incboundary
        BG = boundary(grid(lB), D.D)
        cnt=cnt+1
        ops[cnt] = grid_evaluation_operator(D.S,DiscreteGridSpace(BG,eltype(D.S)),BG)
    end
    for i = 1:length(D.BCs)
        Ac = evaluation_operator(D.S,D.BCs[i].DG)*(D.BCs[i].diff)*ADiff
        ops[i+cnt]=Ac
    end
    BlockOperator(ops)
end

function rhs(D::DiffEquation; incboundary = false, options...)
    op = operator(D; incboundary=incboundary, options...)
    rhs = Array(Array{eltype(src(op)),1},0)
    G, lB = oversampled_grid(D.D,D.S,D.sampling_factor)
    
    op = grid_evaluation_operator(D.S,DiscreteGridSpace(G,eltype(D.S)),G)
    push!(rhs,sample(G,D.DRhs, eltype(src(op))))
    if incboundary
        BG = boundary(grid(lB), D.D)
        push!(rhs,sample(BG,D.DRhs, eltype(src(op))))
    end
    for BC in D.BCs
        push!(rhs,sample(BC.DG,BC.dRhs, eltype(src(op))))
    end
    MultiArray(rhs)
end

function solve(D::DiffEquation, solver=FE_ProjectionSolver; options...)
    G, lB = oversampled_grid(D.D,D.S,D.sampling_factor)
    Adiff= inv_diagonal(D.Diff)
    b = rhs(D; options...)
    OP = operator(D; options...)
    A = solver(OP, length(lB); options...)
    coef  = A * b
    SetFun(D.D, dest(A), Adiff*coef)
end
