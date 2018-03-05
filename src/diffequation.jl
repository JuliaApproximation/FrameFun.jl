# diffequation.jl

"""
A DiffEquation describes a differential equation, with or without boundary conditions.
Parameters:
- Fun is the DictFun that will describe the result

When the equation is solved the equations:
-Diff(Fun) = DRhs on the domain
-diff(Fun) = dRhs on boundary
will hold.
"""

struct BoundaryCondition
    S      :: Dictionary
    diff   :: AbstractOperator
    DG     :: AbstractGrid
end

struct DirichletBC
    dRhs   :: Function
    D      :: Domain
    factor :: Number
    function DirichletBC(dRhs=default_boundary_condition :: Function, D=FullSpace(), factor=1.0)
        new(dRhs,D,factor)
    end
end

struct NeumannBC
    dRhs   :: Function
    D      :: Domain
    function NeumannBC(dRhs=default_boundary_condition :: Function, D=FullSpace())
        new(dRhs,D)
    end
end

BoundaryCondition(S :: Span, D::Domain) = BoundaryCondition(S,IdentityOperator(S),boundary(grid(S),D),default_boundary_condition)
BoundaryCondition(S :: Span, diff::AbstractOperator, D::Domain) = BoundaryCondition(S,diff,boundary(grid(S),D),default_boundary_condition)
BoundaryCondition(S :: Span, diff::AbstractOperator, D::Domain, dRhs::Function) = BoundaryCondition(S,diff,boundary(grid(S),D),dRhs)
default_boundary_condition(x) = 0
default_boundary_condition(x,y) = 0
default_boundary_condition(x,y,z) = 0

function operator(BC :: DirichletBC, S::Span, G::AbstractGrid, D::Domain)
    G = subgrid(G,BC.D)
    BC.factor*grid_evaluation_operator(S,gridspace(G,coeftype(S)),G)
end

function operator(BC :: NeumannBC, S::Span2d, G::AbstractGrid, D::Domain2d)
    G = subgrid(G,BC.D)
    GE = grid_evaluation_operator(S,gridspace(G,coeftype(S)),G)
    dx = Float64[]
    dy = Float64[]
    for i=1:length(G)
        push!(dx, normal(G[i],D)[1])
        push!(dy, normal(G[i],D)[2])
    end
    X = DiagonalOperator(dest(GE), dx)*GE*differentiation_operator(S,(1,0))
    Y = DiagonalOperator(dest(GE), dy)*GE*differentiation_operator(S,(0,1))
    X + Y
end

function operator(BC :: NeumannBC, S::Span1d, G::AbstractGrid1d, D::Domain1d)
    G = subgrid(G,BC.D)
    GE = grid_evaluation_operator(S,gridspace(G,coeftype(S)),G)
    dx = Float64[]
    for i=1:length(G)
        push!(dx, normal(G[i],D)[1])
    end
    DiagonalOperator(dest(GE), dx)*GE*differentiation_operator(S,1)
end

struct DiffEquation
    S     :: Span
    D     :: Domain
    Diff  :: AbstractOperator
    DRhs   :: Function
    BCs    :: Tuple
    sampling_factor
    function DiffEquation(S::Span, D::Domain,Diff::AbstractOperator, DRhs:: Function, BCs::Tuple, sampling_factor=2)
        new(S,D,Diff,DRhs,BCs, sampling_factor)
    end
end

DiffEquation(S::Span, D::Domain, Diff::AbstractOperator, DRhs::Function, BC::BoundaryCondition, sampling_factor=2) = DiffEquation(S,D,Diff,DRhs,(BC,), sampling_factor)

function boundarygrid(D::DiffEquation)
    G, lB = oversampled_grid(D.D,dictionary(D.S),D.sampling_factor)
    boundary(grid(lB),D.D)
end


function operator(D::DiffEquation; incboundary=false, options...)
    #problem = FE_DiscreteProblem(D.D,D.S,2)
    B = D.S
    ADiff = inv_diagonal(D.Diff)
    ops = incboundary ? Array{AbstractOperator}(length(D.BCs)+2,1) : Array{AbstractOperator}(length(D.BCs)+1,1)
    G, lB = oversampled_grid(D.D,dictionary(D.S),D.sampling_factor)

    op = grid_evaluation_operator(D.S,gridspace(G,coeftype(D.S)),G)

    cnt=1
    ops[cnt] = op*D.Diff*ADiff
    BG = boundarygrid(D)
    if incboundary
        cnt=cnt+1
        ops[cnt] = grid_evaluation_operator(D.S,gridspace(BG,coeftype(D.S)),BG)
    end
    for i = 1:length(D.BCs)
        Ac = operator(D.BCs[i],D.S,BG,D.D)*ADiff
        ops[i+cnt]=Ac
    end
    BlockOperator(ops)
end

function rhs(D::DiffEquation; incboundary = false, options...)
    op = operator(D; incboundary=incboundary, options...)
    rhs = Array{Array{coeftype(src(op)),1}}(0)
    G, lB = oversampled_grid(D.D,dictionary(D.S),D.sampling_factor)

    op = grid_evaluation_operator(D.S,gridspace(G,coeftype(D.S)),G)
    push!(rhs,sample(G,D.DRhs, coeftype(src(op))))
    BG = boundarygrid(D)

    if incboundary
        push!(rhs,sample(BG,D.DRhs, coeftype(src(op))))
    end
    for i = 1:length(D.BCs)
        push!(rhs,sample(subgrid(BG,D.BCs[i].D),D.BCs[i].dRhs, coeftype(src(op))))
    end
    MultiArray(rhs)
end

function solve(D::DiffEquation; solver=FE_ProjectionSolver, options...)
    G, lB = oversampled_grid(D.D,dictionary(D.S),D.sampling_factor)
    Adiff= inv_diagonal(D.Diff)
    b = rhs(D; options...)
    OP = operator(D; options...)
    A = solver(OP; scaling=length(lB), options...)
    coef  = A * b
    DictFun(D.D, dictionary(dest(A)), Adiff*coef)
end
