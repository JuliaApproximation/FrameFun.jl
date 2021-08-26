module DiffEquations

using BasisFunctions, GridArrays, DomainSets, LinearAlgebra
using FrameFun
import FrameFun: AZ_A, AZ_Zt

import BasisFunctions: operator, coefficienttype, src, dest, apply!, grid

using BasisFunctions: evaluation

"""
A DiffEquation describes a differential equation, with or without boundary conditions.
Parameters:
- Fun is the Expansion that will describe the result

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

export DirichletBC
struct DirichletBC
    dRhs   :: Function
    D      :: Domain
    factor :: Number
    function DirichletBC(dRhs=default_boundary_condition :: Function, D=DomainSets.FullSpace{Float64}(), factor=1.0)
        new(dRhs,D,factor)
    end
end

export NeumannBC
struct NeumannBC
    dRhs   :: Function
    D      :: Domain
    function NeumannBC(dRhs=default_boundary_condition :: Function, D=DomainSets.FullSpace{Float64}())
        new(dRhs,D)
    end
end

function operator(BC :: DirichletBC, S::Dictionary, G::AbstractGrid, D::Domain)
    G = subgrid(G,BC.D)
    BC.factor*evaluation(coefficienttype(S), S, G)
end

function operator(BC :: NeumannBC, S::Dictionary2d, G::AbstractGrid, D::Domain2d)
    G = subgrid(G,BC.D)
    GE = evaluation(coefficienttype(S), S, G)
    dx = Float64[]
    dy = Float64[]
    for i=1:length(G)
        push!(dx, normal(D,G[i])[1])
        push!(dy, normal(D,G[i])[2])
    end
    X = DiagonalOperator(dest(GE), dx)*GE*differentiation(S,(1,0))
    Y = DiagonalOperator(dest(GE), dy)*GE*differentiation(S,(0,1))
    X + Y
end

function operator(BC :: NeumannBC, S::Dictionary1d, G::AbstractGrid1d, D::Domain1d)
    G = subgrid(G,BC.D)
    GE = evaluation(coefficienttype(S), S, G)
    dx = Float64[]
    for i=1:length(G)
        push!(dx, normal(D,G[i])[1])
    end
    DiagonalOperator(dest(GE), dx)*GE*differentiation(S,1)
end

export DiffEquation
struct DiffEquation
    S     :: Dictionary
    D     :: Domain
    Diff  :: DictionaryOperator
    DRhs   :: Function
    BCs    :: Tuple
    SMP     ::  AbstractOperator
    oversamplingfactor
    function DiffEquation(S::Dictionary, D::Domain,Diff::DictionaryOperator, DRhs:: Function, BCs::Tuple, oversamplingfactor=2)
        SMP = samplingoperator(S, D, oversamplingfactor=oversamplingfactor)
        new(S, D, Diff, DRhs, BCs, SMP, oversamplingfactor)
    end
end

# DiffEquation(S::Dictionary, D::Domain, Diff::DictionaryOperator, DRhs::Function, BC::BoundaryCondition, oversamplingfactor=2) = DiffEquation(S,D,Diff,DRhs,(BC,), oversamplingfactor)


boundarygrid(D::DiffEquation) = boundary(supergrid(grid(dest(D.SMP))),D.D)

function operator(D::DiffEquation; incboundary=false, options...)
    B = D.S
    ADiff = pinv(D.Diff)
    ops = incboundary ? Array{DictionaryOperator}(undef, length(D.BCs)+2,1) : Array{DictionaryOperator}(undef, length(D.BCs)+1,1)
    G = grid(dest(D.SMP))

    op = evaluation(coefficienttype(D.S), D.S, G)
    if length(size(dest(op))) > 1
        op = BasisFunctions.LinearizationOperator(dest(op))*op
    end

    cnt=1
    ops[cnt] = op*D.Diff*ADiff
    BG = boundarygrid(D)
    if incboundary
        cnt=cnt+1
        ops[cnt] = evaluation(coefficienttype(D.S), D.S, BG)
    end
    for i = 1:length(D.BCs)
        Ac = operator(D.BCs[i],D.S,BG,D.D)*ADiff
        ops[i+cnt]=Ac
    end
    BlockOperator(ops)
end

function rhs(D::DiffEquation; incboundary = false, options...)
    op = operator(D; incboundary=incboundary, options...)
    rhs = Array{Array{coefficienttype(src(op)),1}}(undef,0)
    G = grid(dest(D.SMP))

    op = evaluation(coefficienttype(D.S), D.S, G)
    S1 = sample(G,D.DRhs, coefficienttype(src(op)))
    push!(rhs,reshape(S1,length(S1)))
    BG = boundarygrid(D)

    if incboundary
        push!(rhs,sample(BG,D.DRhs, coefficienttype(src(op))))
    end
    for i = 1:length(D.BCs)
        push!(rhs,sample(subgrid(BG,D.BCs[i].D),D.BCs[i].dRhs, coefficienttype(src(op))))
    end
    BlockVector(rhs...)
end

struct PDEApproximation <: ApproximationProblem
    D   ::  DiffEquation
end
coefficienttype(ap::PDEApproximation) = coefficienttype(ap.D.S)

AZ_A(pstyle::ProblemStyle, ap::PDEApproximation; options...) =
        _AZ_A(pstyle, ap; options...)

function AZ_Zt(pstyle::ProblemStyle, ap::PDEApproximation; options...)
    D = ap.D
    G = grid(dest(D.SMP))
    OP = operator(D; options...)
    Zt = 1/length(supergrid(G))*OP'
    Zt
end

import FrameFun: solve
function solve(D::DiffEquation; solverstyle=AZStyle(), options...)
    Adiff = pinv(D.Diff)
    b = rhs(D; options...)
    OP = operator(D; options...)
    A = solver(solverstyle, PDEApproximation(D), OP; options...)
    coef  = A * b
    Expansion(extensionframe(D.D, dest(A)), Adiff*coef)
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


grid_multiplication_operator(a::Function, gb::GridBasis, T=BasisFunctions.operatoreltype(GB)) =
	DiagonalOperator{T}(GB, GB, map(a,grid(GB)))
grid_multiplication_opearator(a::Function, grid::AbstractGrid; options...) =
	grid_multiplication_operator(a,GridBasis(grid); options...)

## Constructors
# All input format checks are done in the first constructor
FECollocationOperator(feframe::ExtensionFrame,pD::Vector,aX::Vector,sampler::DictionaryOperator) = FECollocationOperator{eltype(sampler)}(feframe,pD,aX,sampler)
function FECollocationOperator(feframe::ExtensionFrame,pd::Vector{S1},ax::Vector{S2},samplingfactor::Real) where {S1<:Function,S2<:Function}
    pD = map(p->pseudodifferential_operator(feframe,p),pd)
    gridbasis = GridBasis{coefficienttype(feframe)}(oversampledgrid(feframe,round(Int,samplingfactor*length(feframe))))
    aX = map(a->grid_multiplication_operator(a,gridbasis),ax)
    sampler = evaluation(feframe,grid(gridbasis))
    FECollocationOperator(feframe,pD,aX,sampler)
end
FECollocationOperator(feframe::ExtensionFrame,pd::Vector{S1},ax::Vector{S2}) where {S1<:Function,S2<:Function} = FECollocationOperator(feframe,pd,ax,2)
FECollocationOperator(dom::Domain,ambientdom::Domain,n::Int,pd::Vector{S1},ax::Vector{S2},samplingfactor::Real) where {S1<:Function,S2<:Function} = FECollocationOperator(ExtensionFrame(dom,Fourier(n) → ambientdom),pd,ax,samplingfactor)
FECollocationOperator(dom::Domain,ambientdom::Domain,n,pd::Vector{S1},ax::Vector{S2}) where {S1<:Function,S2<:Function} = FECollocationOperator(dom,ambientdom,n,pd,ax,2)

# Constructors for higher dimensional domains:
FECollocationOperator(dom::Domain,ambientdom::ProductDomain,nn::Tuple,pd::Vector{S1},ax::Vector{S2},samplingfactor::Real) where {S1<:Function,S2<:Function} = FECollocationOperator(ExtensionFrame(dom,tensorproduct(map((x,y)->Fourier(x) → y,nn, components(ambientdom)) )),pd,ax,samplingfactor)
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

end # module
