using BasisFunctions
using StaticArrays
using FrameFun
using Domains



function approximation_step_operators(domain, N1, N2, degr1, degr2)
    B1 = BSplineTranslatesBasis(N1,degr1);
    g1 = BasisFunctions.oversampled_grid(B1,2)
    B2 = BSplineTranslatesBasis(N2,degr2);
    g2 = BasisFunctions.oversampled_grid(B2,2)
    B = B1⊗B2;
    g = g1×g2;

    E1 = CirculantOperator(evaluation_matrix(B1[1],g1)[:])*IndexExtensionOperator(span(B1),gridspace(g1),1:2:length(g1))
    E2 = CirculantOperator(evaluation_matrix(B2[1],g2)[:])*IndexExtensionOperator(span(B2),gridspace(g2),1:2:length(g2))

    G1 = CirculantOperator(E1'E1*[1,zeros(length(g1)-1)...]);
    G2 = CirculantOperator(E2'E2*[1,zeros(length(g2)-1)...]);
    G = G1⊗G2
    E = E1⊗E2

    omega_grid = FrameFun.subgrid(g, domain);
    S = evaluation_operator(span(B), omega_grid).operators[end]

    delta_omega_grid = boundary_grid(g, domain);
    # SO2B = restriction_operator(omega_grid, delta_omega_grid);

    boundary_support = boundary_support_grid(B, delta_omega_grid, omega_grid);
    SO2S = restriction_operator(omega_grid, boundary_support);

    BEO = boundary_extension_operator(delta_omega_grid, B);

    DG = inv(G);
    Eo = S*E
    DEo = Eo*DG

    A = Eo;
    Z = DEo;

    I = IdentityOperator(dest(A));
    P = I-A*Z';
    omega_grid, A, Z, P, SO2S, BEO
end

function system_solve(A, b, epsilon)
    reshape(LAPACK.gelsy!(matrix(A),b,epsilon)[1],size(src(A))...)
end

function projection_solve(A, b, epsilon, R, P, Z)
    W = randn(length(src(A)),R+5);
    AA = apply_multiple(P*A,W)
    y, rr = LAPACK.gelsy!(AA,P*b,epsilon)
    x1 = W*y;
    x1 = reshape(x1, size(src(A))...)
    x2 = Z'*(b-A*x1);
    x1+x2;
end

function selection_solve(A, b, epsilon, P, Z, E, R)
    AA = matrix(R*P*A*E)
    y1, rr = LAPACK.gelsy!(AA,R*P*b,epsilon)
    AA = matrix(R*P*A*E);println(norm(AA*y1-R*P*b))
    x1 = E*y1
    x2 = Z'*(b-A*x1);
    x1+x2;
end

f = (x,y) -> x*(y-1)^2
center = @SVector [.5,.5]
domain  = disk(.3,center)

epsilon=1e-10
N1 = 40; N2 = 40;
degr1 = 1; degr2 = 2;

omega_grid, A, Z, P, SO2S, BEO = approximation_step_operators(domain,N1,N2,degr1,degr2)

b = sample(omega_grid, f);
r = sum(svd(matrix(P*A))[2] .> epsilon)


system_solve(A, b, epsilon)

projection_solve(A, b, epsilon, r, P, Z)

selection_solve(A, b, epsilon, P, Z, BEO, SO2S)
