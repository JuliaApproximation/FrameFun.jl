
# fastsolver.jl

"""
A Randomized Truncated SVD solver storing a decomposition of an operator A.
If A is MxN, and of rank R, the cost of constructing this solver is MR^2.

For tensor product operators it returns a decomposition of the linearized system
"""
immutable TruncatedSvdSolver{ELT} <: AbstractOperator{ELT}
    # Keep the original operator
    op          ::  AbstractOperator
    # Random matrix
    W           ::  MultiplicationOperator
    # Decomposition
    Ut          ::  Array{ELT,2}
    VS          ::  Array{ELT,2}
    # For storing intermediate results when applying
    y           ::  Array{ELT,1}
    sy          ::  Array{ELT,1}
    function TruncatedSvdSolver(op::AbstractOperator; cutoff = default_cutoff(problem), rank_estimate = 5, verbose = false, options...)
        finished=false
        USV = ()
        R = rank_estimate
        while R<=size(op,2)
            random_matrix = map(ELT, rand(size(op,2), R))
            Wsrc = ELT <: Complex ? Cn{ELT}(size(random_matrix,2)) : Rn{ELT}(size(random_matrix,2))
            Wdest = src(op)
            W = MatrixOperator(Wsrc, Wdest, random_matrix)
        
            USV = LAPACK.gesdd!('S',matrix(op * W))
            verbose && println("minimal singular value: ",minimum(USV[2])," cutoff : ",cutoff)
            if minimum(USV[2])<cutoff || R==size(op,2)
                S = USV[2]
                maxind = findlast(S.>cutoff)
                Sinv = 1./S[1:maxind]
                y = zeros(ELT, size(USV[3],1))
                sy = zeros(ELT, maxind)
                return new(op, W, USV[1][:,1:maxind]',USV[3][1:maxind,:]'*diagm(Sinv[:]),y,sy)
            else
                R = min(round(Int,2*R),size(op,2))
            end
        end
    end
end

src(t::TruncatedSvdSolver) = dest(t.op)
dest(t::TruncatedSvdSolver) = src(t.op)
inv(t::TruncatedSvdSolver) = t.op

TruncatedSvdSolver(op::AbstractOperator; options...) =
    TruncatedSvdSolver{eltype(op)}(op::AbstractOperator; options...)

function apply!(s::TruncatedSvdSolver, destset, srcset, coef_dest, coef_src)
    # Applying the truncated SVD to P*Rhs
    A_mul_B!(s.sy, s.Ut, coef_src)
    A_mul_B!(s.y, s.VS, s.sy)
    apply!(s.W, coef_dest, s.y)
end
