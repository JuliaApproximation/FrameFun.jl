
function FNA2_kappa_lambda(platform, N, M, threshold; options...)
    local G_N

    # Computation of the Gram matrix may be expensive: check if it was passed as an argument
    if haskey(options, :G_N)
        G_N = options[:G_N]
    else
        G_N = matrix(discretization(platform, N, samplingstyle=GramStyle(); options...))
    end

    G_MN = solver(platform, N, samplingstyle=OversamplingStyle(), M=M, solverstyle=DirectStyle(), directsolver=:svd, warnslow=false, options...)
    U = G_MN.solver.U
    S = G_MN.solver.S
    V = collect(G_MN.solver.V)
    I = findlast(S .>= threshold)
    U_reg = U[:,1:I]
    V_reg = V[:,1:I]
    S_reg = S[1:I]

    # Computation of kappa
    B_MN = V_reg * inv(Diagonal(S_reg)) * U_reg'
    Q = svdvals(B_MN' * G_N * B_MN)
    kappa = sqrt(Q[1])

    # Computation of lambda
    V_eps = V[:,I+1:end]
    C_MN = V_eps * V_eps'
    Q2 = svdvals(C_MN' * G_N * C_MN)
    lambda = sqrt(Q2[1])/threshold

    kappa, lambda
end
