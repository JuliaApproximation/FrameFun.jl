
"""
A fast FE solver based on a low-rank approximation of the plunge region. The plunge region
is isolated using a projection operator. This algorithm contains an extra smoothing step
"""
struct AZSmoothSolver{T} <: AbstractSolverOperator{T}
    TS          ::  DictionaryOperator
    A           ::  DictionaryOperator
    Zt          ::  DictionaryOperator
    plunge_op   ::  DictionaryOperator    # store the operator because it allocates memory
    b           ::  Array{T,1}
    blinear     ::  Array{T,1}
    syv         ::  Array{T,1}
    x1          ::  Array{T}
    x2          ::  Array{T}
    x3          ::  Array{T}
    Q           ::  Array{T,2}
    D           ::  DictionaryOperator
    AD          ::  DictionaryOperator

    function AZSmoothSolver{T}(A, Zt, D, plunge_op, AD, TS1, TS2, Q) where {T}
        # Pre-allocation
        b = zeros(size(dest(plunge_op)))
        blinear = zeros(T, length(dest(A)))
        x1 = zeros(size(src(A)))
        x2 = zeros(size(src(A)))
        x3 = zeros(size(src(A)))
        syv = zeros(size(TS2.Ut,1))
        new(TS1, A, Zt, plunge_op, b, blinear, syv, x1, x2, x3, Q, D, AD)
    end
end

function AZSmoothSolver(A::DictionaryOperator, Zt::DictionaryOperator; options...)
    D = WeightedSmoothingOperator(src(A); options...)
    AZSmoothSolver(A, Zt, D; options...)
end

function AZSmoothSolver(A::DictionaryOperator{T}, Zt::DictionaryOperator, D::DictionaryOperator;
            REG = default_regularization,
            rankestimate = 40,
            threshold = default_threshold(A),
            thresholdv = sqrt(threshold),
            options...) where {T}

    plunge_op = plunge_operator(A, Zt)
    TS1 = REG(plunge_op*A; threshold = threshold, rankestimate=rankestimate, options...)
    TS2 = REG(Zt*plunge_op; threshold = thresholdv, rankestimate=rankestimate, options...)
    AD = inv(D)
    ADV = (TS2.Ut)'.*diagonal(AD)
    # Orthogonal basis for D^(-1)V_mid
    Q, R = qr(ADV)
    Q = Matrix(Q)
    AZSmoothSolver{T}(A, Zt, D, plunge_op, AD, TS1, TS2, Q)
end


operator(s::AZSmoothSolver) = s.A

apply!(s::AZSmoothSolver, coef_dest, coef_src) =
    _apply!(s, coef_dest, coef_src, s.TS, s.A, s.Zt, s.plunge_op, s.b, s.blinear, s.syv, s.x1, s.x2, s.x3, s.Q, s.D, s.AD)

function _apply!(s::AZSmoothSolver, coef_dest, coef_src, TS, A, Zt, plunge_op, b, blinear, syv, x1, x2, x3, Q, D, AD)
    # Apply plunge operator to right hand side
    apply!(plunge_op, b, coef_src)
    BasisFunctions.linearize_coefficients!(dest(A), blinear, b)
    # x2 is solving for the middle singular values
    apply!(TS, x2, blinear)
    # Project Dx2 onto the orthogonal complement of D^(-1)V_mid
    apply!(D, x1, x2)
    ## apply!(MantrixOperator(s.Q'),s.syv,s.x1)
    ## apply!(MatrixOperator(s.Q),s.x3,s.syv)
    ## for i = 1:length(s.x1)
    ##     s.x1[i] = s.x1[i]- s.x3[i]
    ## end
    ## apply!(AD,s.x3,s.x1)
    ## for i = 1:length(s.x2)
    ##     s.x2[i] = s.x2[i] - s.x3[i]
## end
    apply!(MatrixOperator(Matrix(Q')), syv, x1)
    apply!(MatrixOperator(Q), x3, syv)
    apply!(AD, x2, x3)

    # post smoothing step
    apply!(A, b, x2)
    apply!(Zt, x1, coef_src-b)
    for i in 1:length(coef_dest)
        coef_dest[i] = x1[i] + x2[i]
    end
    coef_dest
end

# For Dictionary's that have a DC component
dc_index(b::ChebyshevBasis) = 1
dc_index(b::FourierBasis) = 1

function WeightedSmoothingOperator(dict::Dictionary; scaling = default_scaling_function, order = 1, options...)
    if order == 1
        WeightedSmoothingOperator(dict, scaling)
    else
        WeightedSmoothingOperator(dict, (d,i) -> scaling(d,i)^order)
    end
end

WeightedSmoothingOperator(dict::Dictionary, scaling) = DiagonalOperator(dict, dict, [scaling(dict, i) for i in eachindex(dict)][:])

function default_scaling_function(dict::Dictionary1d, idx)
    f = abs(value(native_index(dict, idx)))
    1 + f
end

default_scaling_function(dict::TensorProductDict, I) =
    default_scaling_function(element(dict, 1), I[1]) + default_scaling_function(element(dict, 2), I[2])
