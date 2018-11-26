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
    diff   :: DictionaryOperator
    DG     :: AbstractGrid
end

struct DirichletBC
    dRhs   :: Function
    D      :: Domain
    factor :: Number
    function DirichletBC(dRhs=default_boundary_condition :: Function, D=DomainSets.FullSpace{Float64}(), factor=1.0)
        new(dRhs,D,factor)
    end
end

struct NeumannBC
    dRhs   :: Function
    D      :: Domain
    function NeumannBC(dRhs=default_boundary_condition :: Function, D=DomainSets.FullSpace{Float64}())
        new(dRhs,D)
    end
end

# BoundaryCondition(S :: Dictionary, D::Domain) = BoundaryCondition(S,IdentityOperator(S),boundary(grid(S),D),default_boundary_condition)
# BoundaryCondition(S :: Dictionary, diff::DictionaryOperator, D::Domain) = BoundaryCondition(S,diff,boundary(grid(S),D),default_boundary_condition)
# BoundaryCondition(S :: Dictionary, diff::DictionaryOperator, D::Domain, dRhs::Function) = BoundaryCondition(S,diff,boundary(grid(S),D),dRhs)
# default_boundary_condition(x) = 0
# default_boundary_condition(x,y) = 0
# default_boundary_condition(x,y,z) = 0

function operator(BC :: DirichletBC, S::Dictionary, G::AbstractGrid, D::Domain)
    G = subgrid(G,BC.D)
    BC.factor*grid_evaluation_operator(S,gridbasis(G,coeftype(S)),G)
end

function operator(BC :: NeumannBC, S::Dictionary2d, G::AbstractGrid, D::Domain2d)
    G = subgrid(G,BC.D)
    GE = grid_evaluation_operator(S,gridbasis(G,coeftype(S)),G)
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

function operator(BC :: NeumannBC, S::Dictionary1d, G::AbstractGrid1d, D::Domain1d)
    G = subgrid(G,BC.D)
    GE = grid_evaluation_operator(S,gridbasis(G,coeftype(S)),G)
    dx = Float64[]
    for i=1:length(G)
        push!(dx, normal(G[i],D)[1])
    end
    DiagonalOperator(dest(GE), dx)*GE*differentiation_operator(S,1)
end

struct DiffEquation
    S     :: Dictionary
    D     :: Domain
    Diff  :: DictionaryOperator
    DRhs   :: Function
    BCs    :: Tuple
    sampling_factor
    function DiffEquation(S::Dictionary, D::Domain,Diff::DictionaryOperator, DRhs:: Function, BCs::Tuple, sampling_factor=2)
        new(S,D,Diff,DRhs,BCs, sampling_factor)
    end
end

# DiffEquation(S::Dictionary, D::Domain, Diff::DictionaryOperator, DRhs::Function, BC::BoundaryCondition, sampling_factor=2) = DiffEquation(S,D,Diff,DRhs,(BC,), sampling_factor)

function boundarygrid(D::DiffEquation)
    G, lB = oversampled_grid(D.D,D.S,D.sampling_factor)
    boundary(grid(lB),D.D)
end


function operator(D::DiffEquation; incboundary=false, options...)
    B = D.S
    ADiff = pinv(D.Diff)
    ops = incboundary ? Array{DictionaryOperator}(undef, length(D.BCs)+2,1) : Array{DictionaryOperator}(undef, length(D.BCs)+1,1)
    G, lB = oversampled_grid(D.D,D.S,D.sampling_factor)

    op = grid_evaluation_operator(D.S,gridbasis(G,coeftype(D.S)),G)

    cnt=1
    ops[cnt] = op*D.Diff*ADiff
    BG = boundarygrid(D)
    if incboundary
        cnt=cnt+1
        ops[cnt] = grid_evaluation_operator(D.S,gridbasis(BG,coeftype(D.S)),BG)
    end
    for i = 1:length(D.BCs)
        Ac = operator(D.BCs[i],D.S,BG,D.D)*ADiff
        ops[i+cnt]=Ac
    end
    BlockOperator(ops)
end

function rhs(D::DiffEquation; incboundary = false, options...)
    op = operator(D; incboundary=incboundary, options...)
    rhs = Array{Array{coeftype(src(op)),1}}(undef,0)
    G, lB = oversampled_grid(D.D,D.S,D.sampling_factor)

    op = grid_evaluation_operator(D.S,gridbasis(G,coeftype(D.S)),G)
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

function solve(D::DiffEquation; solver=AZSolver, options...)
    G, lB = oversampled_grid(D.D,D.S,D.sampling_factor)
    Adiff= pinv(D.Diff)
    b = rhs(D; options...)
    OP = operator(D; options...)
    A = solver(OP; scaling=length(lB), options...)
    coef  = A * b
    DictFun(D.D, dest(A), Adiff*coef)
end



""" Spectral collocation Operators
A spectral collocation operator is a DictionaryOperator from a Fourier Extension Frame to a GridBasis.
It represents the linear differential operator:
L[u](x) = aX[1](X)*sampler*pD[1](D) + ⋯ + aX[r](X)*sampler*pD[r](D)
where X is a grid multiplication operator for x -> x, and D is the derivative operator
"""
struct FECollocationOperator{T} <: DictionaryOperator{T}
  # L = aX[1]*sampler*pD[1] + ⋯ + aX[r]*sampler*pD[r]
  feframe :: ExtensionFrame # This should be Fourier Extension of odd order to work
  pD :: Vector
  aX :: Vector
  sampler :: DictionaryOperator{T}
  function FECollocationOperator{T}(feframe::ExtensionFrame,pD::Vector,aX::Vector,sampler::DictionaryOperator{T}) where {T}
      r = length(pD)
      @assert r == length(aX)
      for k = 1:r
          @assert src(pD[k]) == feframe
          @assert dest(pD[k]) == feframe # Assumed feframe is of odd order
          @assert src(sampler) == feframe
          @assert src(aX[k]) == dest(sampler)
      end
      new(feframe,pD,aX,sampler)
  end
end

## Constructors
# All input format checks are done in the first constructor
FECollocationOperator(feframe::ExtensionFrame,pD::Vector,aX::Vector,sampler::DictionaryOperator) = FECollocationOperator{eltype(sampler)}(feframe,pD,aX,sampler)
function FECollocationOperator(feframe::ExtensionFrame,pd::Vector{S1},ax::Vector{S2},samplingfactor::Real) where {S1<:Function,S2<:Function}
    pD = map(p->pseudodifferential_operator(feframe,p),pd)
    gridbasis = GridBasis(oversampled_grid(feframe,samplingfactor)[1],coefficient_type(feframe))
    aX = map(a->grid_multiplication_operator(a,gridbasis),ax)
    sampler = evaluation_operator(feframe,grid(gridbasis))
    FECollocationOperator(feframe,pD,aX,sampler)
end
FECollocationOperator(feframe::ExtensionFrame,pd::Vector{S1},ax::Vector{S2}) where {S1<:Function,S2<:Function} = FECollocationOperator(feframe,pd,ax,2)
FECollocationOperator(dom::Domain,ambientdom::Domain,n::Int,pd::Vector{S1},ax::Vector{S2},samplingfactor::Real) where {S1<:Function,S2<:Function} = FECollocationOperator(ExtensionFrame(dom,FourierBasis(n,leftendpoint(ambientdom),rightendpoint(ambientdom))),pd,ax,samplingfactor)
FECollocationOperator(dom::Domain,ambientdom::Domain,n,pd::Vector{S1},ax::Vector{S2}) where {S1<:Function,S2<:Function} = FECollocationOperator(dom,ambientdom,n,pd,ax,2)

# Constructors for higher dimensional domains:
FECollocationOperator(dom::Domain,ambientdom::ProductDomain,nn::Tuple,pd::Vector{S1},ax::Vector{S2},samplingfactor::Real) where {S1<:Function,S2<:Function} = FECollocationOperator(ExtensionFrame(dom,tensorproduct(map((x,y)->FourierBasis(x,leftendpoint(y),rightendpoint(y)),nn, elements(ambientdom)))),pd,ax,samplingfactor)
FECollocationOperator(dom::Domain,ambientdom::ProductDomain,nn::Tuple,pd::Vector{S1},ax::Vector{S2}) where {S1<:Function,S2<:Function} = FECollocationOperator(dom,ambientdom,nn,pd,ax,2)

src(L::FECollocationOperator) = L.feframe
dest(L::FECollocationOperator) = dest(sampler(L))

function apply!(L::FECollocationOperator,dest,src,coef_dest,coef_src)
    coef_dest[:] = L.aX[1]*(L.sampler*(L.pD[1]*coef_src))
    for k = 2:length(L.pD)
        coef_dest[:] += L.aX[k]*(L.sampler*(L.pD[k]*coef_src))
    end
    coef_dest
end

sampler(L) = L.sampler
grid(L::FECollocationOperator) = grid(sampler(L))
